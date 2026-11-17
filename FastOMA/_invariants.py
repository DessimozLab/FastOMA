"""Checks that no protein is lost or duplicated on its way through the pipeline.

A violated invariant is a bug. In production the caller logs it and repairs what it can (e.g. keeps the
protein as a singleton), so a run still finishes. With the environment variable FASTOMA_STRICT_INVARIANTS
set (as the tests do) it raises instead, so the repair cannot hide the bug.
"""
import os


def invariant_violated(log, msg):
    if os.environ.get("FASTOMA_STRICT_INVARIANTS"):
        raise AssertionError(msg)
    log(msg)
