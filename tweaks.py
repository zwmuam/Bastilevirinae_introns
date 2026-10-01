#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "GNU GENERAL PUBLIC LICENSE"
__email__ = "jakub.barylski@gmail.com"

import inspect
import pickle
import sys
import warnings
from functools import wraps
from multiprocessing import cpu_count
from pathlib import Path
from subprocess import run, DEVNULL
from typing import Callable, Dict, Collection, List, Hashable, Any

import joblib
import numpy as np
from loguru import logger
from tqdm import tqdm


# LOGGING format and configuration
def loguru_showwarning(message, category, filename, lineno, file=None, line=None):
    filename = Path(filename).name
    logger.warning(f"{filename} {category.__name__}: {message}")


logger.remove()
log_format = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> <cyan>{function}</cyan>: <level>{message}</level>"
logger.add(sys.stderr, format=log_format)
warnings.showwarning = loguru_showwarning

# DEFAULTS
default_threads = max(cpu_count() - 1, 1)
seed_max = 2 ** 32 - 1

extensions = {'gbk': {'.gb', '.gbk'},
              'gff': {'.gff', '.gff3'},
              'fna': {'.fa', '.fas', '.fasta', '.fna'},
              'faa': {'.fa', '.fas', '.fasta', '.faa'}}

for form in extensions:  # add gzip compressed files
    extensions[form].update([f'{e}.gz' for e in extensions[form]])


def frantic_search(dictionary: Dict[Hashable, Any],
                   *possible_keys: Hashable):
    """
    Find first of the several keys that are in a dictionary and return the value
    :param dictionary: dictionary to search
    :param possible_keys: any number of possible keys (preferred first)
    :return: dictionary[first_found_key]
    """
    for key in possible_keys:
        if key in dictionary:
            return dictionary[key]
    missed_keys = ', '.join([str(k) for k in possible_keys])
    raise KeyError(f'Found none of the: {missed_keys}')


def find_files(directory: Path,
               file_type: str,
               descent: bool = False) -> List[Path]:
    """
    Find files of a given type in a folder.
    :param directory: directory to search in
    :param file_type: file type to search for
    :param descent: whether to search in subdirectories (and their subdirectories)
    :return:
    """
    main_file_catalogue = []
    expected_extensions = extensions[file_type]
    detected_files = [f for f in directory.iterdir() if f.suffix in expected_extensions]
    subdirectories = [f for f in directory.iterdir() if f.is_dir()]

    if detected_files:
        main_file_catalogue.extend(detected_files)

    elif subdirectories and descent:
        for sd in subdirectories:
            main_file_catalogue.extend(find_files(sd, file_type, descent))

    return main_file_catalogue


# PARALLELIZATION
import mpire
from mpire import WorkerPool


class Parallel:
    """
    Modern high-performance multiprocessing wrapper around MPIRE WorkerPool.
    Optimized for multi-core Linux workstations, bioinformatic workflows,
    complex Python/Biopython object graphs, dynamic chunking, and seed safety.
    """

    def __init__(self,
                 parallelized_function: Callable,
                 input_collection: Collection = None,
                 random_replicates: int = None,
                 kwargs: Dict = None,
                 n_jobs: int = None,
                 description: str = None,
                 bar: bool = True,
                 bar_color: str = 'cyan',
                 chunk_size: Any = 'auto',
                 shared_objects: Any = None,
                 keep_order: bool = True,
                 **extra_kwargs):

        assert bool(random_replicates) ^ bool(input_collection is not None), \
            'You need to specify EITHER an input collection OR number of random replicates'

        n_jobs = n_jobs if n_jobs is not None else default_threads
        kwargs = {} if not kwargs else kwargs
        description = description if description else parallelized_function.__name__

        # Handle random replicates vs input collection
        if random_replicates:
            ss = np.random.SeedSequence()
            child_seeds = ss.spawn(random_replicates)
            input_collection = [int(s.generate_state(1)[0]) for s in child_seeds]
            description = f'{description} 🎲'

        input_list = list(input_collection) if not isinstance(input_collection, list) else input_collection
        total_items = len(input_list)

        # Dynamic chunk size estimation for variable size bioinformatic processes
        if chunk_size == 'auto':
            if total_items > 0 and n_jobs > 0:
                chunk_size = max(1, total_items // (n_jobs * 4))
            else:
                chunk_size = 1

        worker_kwargs = {'shared_objects': shared_objects} if shared_objects else {}

        def _worker_wrapper(item):
            return parallelized_function(item, **kwargs)

        with WorkerPool(n_jobs=n_jobs, **worker_kwargs) as pool:
            pbar_opts = {'desc': description, 'colour': bar_color} if bar else {}
            map_func = pool.map if keep_order else pool.map_unordered
            self.result = map_func(_worker_wrapper,
                                   input_list,
                                   progress_bar=bar,
                                   progress_bar_options=pbar_opts if bar else None,
                                   chunk_size=chunk_size)

    def print_progress(self):
        pass


def run_external(command: List[str],
                 stdout='suppress',
                 stdin=None):
    """
    Run external (non-python) command
    :param command: list of the expressions
                    that make up the shell command e.g. ['ls', '-lh']
    :param stdout: do not print the log messages
    :param stdin: input for the command
    """
    sanitized_command = [str(c) for c in command]

    logger.info(" ".join(sanitized_command))
    if stdout == 'suppress':
        process = run(sanitized_command, stdout=DEVNULL, stderr=DEVNULL, input=stdin)
    elif stdout == 'capture':
        process = run(sanitized_command, capture_output=True, input=stdin)
        if process.returncode != 0:
            stderr_msg = process.stderr.decode() if process.stderr else ''
            raise ChildProcessError(f'"{" ".join(sanitized_command)}" crashed with:\n{stderr_msg}')
        return process.stdout
    else:
        process = run(sanitized_command)
    if process.returncode != 0:
        stderr_msg = process.stderr.decode() if getattr(process, 'stderr', None) else ''
        raise ChildProcessError(f'"{" ".join(sanitized_command)}" crashed with:\n{stderr_msg}')


def parse_fasta(fasta: Path):
    """
    Simple and relatively fast fasta parser
    used when no complex sequence handling is required
    :param fasta: path to a fasta file
    :return: generator yielding tuples of (identifier, sequence)
    """
    identifier, sequence = None, []
    with fasta.open() as fas:
        for line in fas:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if identifier is not None:
                    yield identifier, ''.join(sequence)
                identifier = line.lstrip('>').split(' ')[0]
                sequence = []
            else:
                sequence.append(line)
    if identifier is not None:
        yield identifier, ''.join(sequence)


def checkpoint(funct: callable):
    """
    Simple serialization decorator
    that saves the function result
    if exacted output file doesn't exist or is empty
    or read it if it is non-empty
    @param funct: function to wrap
    @param pickle_path: a path to an output file
    @param serialization_method: a module used for serialization (either joblib or pickle)
    @return:
    """

    signature = inspect.signature(funct)

    @wraps(funct)
    def save_checkpoint(*args, **kwargs):

        bound_args = signature.bind(*args, **kwargs)
        pickle_path = Path(bound_args.arguments.get('pickle_path',
                                                    signature.parameters['pickle_path'].default))
        if pickle_path:
            try:
                with pickle_path.open('rb') as file_object:
                    result = pickle.load(file_object)
                logger.info(f'\ntemporary file read from: {pickle_path.as_posix()}\n', flush=True)
                return result
            except (FileNotFoundError, IOError, EOFError):
                sys.setrecursionlimit(5000)
                result = funct(*args, **kwargs)
                with pickle_path.open('wb') as out:
                    pickle.dump(result, out)
                logger.info(f'\ntemporary file stored at: {pickle_path.as_posix()}\n', flush=True)
                return result

    return save_checkpoint
