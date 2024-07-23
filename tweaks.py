import inspect
import pickle
import sys
import warnings
from functools import wraps
from itertools import chain
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
    Find first of several keys that is in a dictionary and return value
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
               descent: bool = False):
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
class Parallel(joblib.Parallel):
    """
    The modification of joblib.Parallel
    with a TQDM progress bar
    according to Nth
    (https://stackoverflow.com/questions/37804279/how-can-we-use-tqdm-in-a-parallel-execution-with-joblib)
    """

    def __init__(self,
                 parallelized_function: Callable,
                 input_collection: Collection = None,
                 random_replicates: int = None,
                 kwargs: Dict = None,
                 n_jobs=None,
                 backend=None,
                 description: str = None,
                 verbose=0,
                 timeout=None,
                 pre_dispatch='2 * n_jobs',
                 batch_size='auto',
                 temp_folder=None, max_nbytes='1M', mmap_mode='r',
                 prefer=None,
                 require=None,
                 bar: bool = True):

        if not n_jobs:
            n_jobs = default_threads

        joblib.Parallel.__init__(self,
                                 n_jobs=n_jobs,
                                 backend=backend,
                                 verbose=verbose,
                                 timeout=timeout,
                                 pre_dispatch=pre_dispatch,
                                 batch_size=batch_size,
                                 temp_folder=temp_folder,
                                 max_nbytes=max_nbytes,
                                 mmap_mode=mmap_mode,
                                 prefer=prefer,
                                 require=require)

        assert bool(random_replicates) ^ bool(input_collection), 'You need to specify EITHER an input collection ' \
                                                                 'OR number of random replicates of the function'

        kwargs = {} if not kwargs else kwargs

        description = description if description else parallelized_function.__name__

        if random_replicates:
            input_collection = np.random.randint(0, 2 ** 32 - 1,
                                                 size=random_replicates,
                                                 dtype=np.int64)
            description = f'{description} 🎲'

        jobs = ((joblib.delayed(parallelized_function)(e, **kwargs)) for e in input_collection)

        if bar:
            self._progress = tqdm(total=len(input_collection), file=sys.stdout)
            if description:
                self._progress.set_description(description)
        else:
            self._progress = None

        self.result = list(self.__call__(jobs))

        if self._progress:
            self._progress.close()
            print(flush=True)

    def print_progress(self):
        if self._progress:
            self._progress.n = self.n_completed_tasks
            self._progress: tqdm
            self._progress.refresh()


class BatchParallel(Parallel):
    """ Version of the Parallel used for large numbers of
    computationally un-intensive processes """

    def __init__(self,
                 parallelized_function: Callable,
                 input_collection: Collection = None,
                 random_replicates: int = None,
                 chunk_size: float = 0.1 / default_threads,
                 kwargs: Dict = {},
                 n_jobs=None,
                 backend=None,
                 description: str = None,
                 verbose=0,
                 timeout=None,
                 pre_dispatch='2 * n_jobs',
                 batch_size='auto',
                 temp_folder=None, max_nbytes='1M', mmap_mode='r',
                 prefer=None,
                 require=None,
                 bar: bool = True):

        assert bool(random_replicates) ^ bool(input_collection), 'You need to specify EITHER an input collection ' \
                                                                 'OR number of random replicates of the function'

        if description is None:
            description = parallelized_function.__name__

        if random_replicates:
            input_collection = np.random.randint(0, 2 ** 32 - 1,
                                                 size=random_replicates,
                                                 dtype=np.int64)
            description = f'{description} 🎲'

        def wrapper_function(batch):
            return tuple([parallelized_function(element, **kwargs) for element in batch])

        chunk_n = int(len(input_collection) * chunk_size)

        batches = [input_collection[i * chunk_n:(i + 1) * chunk_n] for i in
                   range((len(input_collection) + chunk_n - 1) // chunk_n)]

        Parallel.__init__(self,
                          parallelized_function=wrapper_function,
                          input_collection=batches,
                          n_jobs=n_jobs,
                          backend=backend,
                          verbose=verbose,
                          timeout=timeout,
                          pre_dispatch=pre_dispatch,
                          batch_size=batch_size,
                          temp_folder=temp_folder,
                          max_nbytes=max_nbytes,
                          mmap_mode=mmap_mode,
                          prefer=prefer,
                          require=require,
                          description=description,
                          bar=bar)

        self.result = chain.from_iterable(self.result)


def run_external(command: List[str],
                 stdout='suppress',
                 stdin=None):
    """
    Run external (non-python) command
    :param command: list of the expressions
                    that make up shell command e.g. ['ls', '-lh']
    :param stdout: do not print the log messages
    :param stdin: input for the command
    """
    sanitized_command = [str(c) for c in command]

    logger.info(" ".join(sanitized_command))
    if stdout == 'suppress':
        process = run(sanitized_command, stdout=DEVNULL, stderr=DEVNULL, input=stdin)
    elif stdout == 'capture':
        process = run(sanitized_command, capture_output=True, input=stdin)
        return process.stdout
    else:
        process = run(sanitized_command)
    if process.returncode or process.stderr:
        raise ChildProcessError(f'"{" ".join(sanitized_command)}" crashed with:\n{process.stderr}')


def parse_fasta(fasta: Path):
    """
    todo
    :param fasta:
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
    yield identifier, ''.join(sequence)


def checkpoint(funct: callable):
    """
    Simple serialization decorator
    that saves function result
    if exacted output file don't exist or is empty
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
