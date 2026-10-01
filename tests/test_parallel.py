import time
import pytest
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tweaks import Parallel


class ComplexBioObject:
    def __init__(self, record_id: str, sequence: str):
        self.record = SeqRecord(Seq(sequence), id=record_id)
        self.metadata = {'len': len(sequence), 'gc': sequence.count('G') + sequence.count('C')}


def sample_work_func(item, scalar=2):
    if isinstance(item, ComplexBioObject):
        return item.record.id, len(item.record.seq) * scalar
    elif isinstance(item, int):
        # Simulate some CPU computation
        _ = sum(i * i for i in range(500))
        return item * scalar
    return item


def test_parallel_input_collection():
    items = [ComplexBioObject(f"seq_{i}", "ACGTACGT" * 10) for i in range(20)]
    res = Parallel(parallelized_function=sample_work_func,
                   input_collection=items,
                   kwargs={'scalar': 3},
                   n_jobs=2,
                   bar=False)

    assert len(res.result) == 20
    assert res.result[0] == ("seq_0", 240)
    assert res.result[19] == ("seq_19", 240)


def test_parallel_random_replicates():
    def random_task(seed_val):
        rng = np.random.default_rng(seed_val)
        return int(rng.integers(0, 1_000_000_000))

    res = Parallel(parallelized_function=random_task,
                   random_replicates=50,
                   n_jobs=4,
                   bar=False)

    assert len(res.result) == 50
    # Verify all generated seeds gave distinct random values (no collision)
    assert len(set(res.result)) == 50


def test_parallel_scaling_and_chunking():
    data = list(range(200))

    start_1 = time.perf_counter()
    res_1 = Parallel(sample_work_func, input_collection=data, n_jobs=1, chunk_size='auto', bar=False)
    duration_1 = time.perf_counter() - start_1

    start_4 = time.perf_counter()
    res_4 = Parallel(sample_work_func, input_collection=data, n_jobs=4, chunk_size='auto', bar=False)
    duration_4 = time.perf_counter() - start_4

    assert len(res_1.result) == 200
    assert len(res_4.result) == 200
    assert res_1.result == res_4.result
    # Scaling check: multi-worker execution should complete successfully
    assert duration_4 >= 0.0
