import pytest
from pathlib import Path
import tempfile
from tweaks import frantic_search, find_files, parse_fasta, run_external


def test_frantic_search():
    d = {'a': 1, 'b': 2, 'c': 3}
    assert frantic_search(d, 'x', 'y', 'b', 'c') == 2
    assert frantic_search(d, 'c', 'a') == 3
    with pytest.raises(KeyError):
        frantic_search(d, 'x', 'z')


def test_find_files(tmp_path):
    f1 = tmp_path / "test1.fasta"
    f2 = tmp_path / "test2.fna"
    f3 = tmp_path / "test3.txt"
    f1.write_text(">seq1\nACGT\n")
    f2.write_text(">seq2\nTGCA\n")
    f3.write_text("hello\n")

    sub_dir = tmp_path / "sub"
    sub_dir.mkdir()
    f4 = sub_dir / "test4.fa"
    f4.write_text(">seq4\nAAAA\n")

    files = find_files(tmp_path, 'fna', descent=False)
    file_names = {f.name for f in files}
    assert "test1.fasta" in file_names or "test2.fna" in file_names

    files_descent = find_files(tmp_path, 'fna', descent=True)
    file_names_descent = {f.name for f in files_descent}
    assert "test1.fasta" in file_names_descent or "test2.fna" in file_names_descent


def test_parse_fasta(tmp_path):
    f = tmp_path / "example.fasta"
    f.write_text(">seq1 first sequence\nACGT\nACGT\n>seq2 second sequence\nGGGG\n")

    records = list(parse_fasta(f))
    assert len(records) == 2
    assert records[0] == ("seq1", "ACGTACGT")
    assert records[1] == ("seq2", "GGGG")

    empty_f = tmp_path / "empty.fasta"
    empty_f.write_text("")
    assert list(parse_fasta(empty_f)) == []


def test_run_external():
    # Simple echo test
    output = run_external(["echo", "hello"], stdout='capture')
    assert output.strip() == b"hello"

    # Command returning error code
    with pytest.raises(ChildProcessError):
        run_external(["ls", "/non_existent_file_path_12345"], stdout='capture')
