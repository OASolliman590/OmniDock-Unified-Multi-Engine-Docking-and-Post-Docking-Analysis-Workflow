import types

import pytest

from workflow import execution


@pytest.mark.parametrize(
    "successful,total,status",
    [(3, 3, "completed"), (1, 3, "partial"), (0, 3, "failed"), (0, 0, "failed")],
)
def test_batch_status_reflects_failures(monkeypatch, tmp_path, successful, total, status):
    fake_module = types.ModuleType("batch_pdb_preparation")

    class FakePipeline:
        def __init__(self, *_args, **_kwargs):
            pass

        def run_batch_processing(self):
            return {"successful": successful, "total": total}

    fake_module.BatchPDBPreparationPipeline = FakePipeline
    monkeypatch.setitem(__import__("sys").modules, "batch_pdb_preparation", fake_module)
    result = execution.run_pdb_batch(str(tmp_path / "batch.yaml"), str(tmp_path / "out"))
    assert result.status == status


def test_cli_exits_non_zero_unless_every_entry_succeeded(monkeypatch, tmp_path):
    from workflow import cli

    for status, code in (("completed", 0), ("partial", 1), ("failed", 1)):
        monkeypatch.setattr(execution, "run_pdb_batch", lambda *_a, _s=status, **_k: types.SimpleNamespace(status=_s))
        assert cli.main(["pdb", "batch", "--config", str(tmp_path / "c.yaml"), "--output", str(tmp_path / "o")]) == code
