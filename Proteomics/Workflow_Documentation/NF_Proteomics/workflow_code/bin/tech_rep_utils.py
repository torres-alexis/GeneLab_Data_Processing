"""Shared helpers for Has Tech Reps / Source Name technical-replicate handling."""

from __future__ import annotations

import re


def parse_has_tech_reps(val) -> bool | None:
    if val is None:
        return None
    s = str(val).strip()
    if not s:
        return None
    u = s.upper()
    if u in ("TRUE", "T", "1", "YES"):
        return True
    if u in ("FALSE", "F", "0", "NO"):
        return False
    return None


def fragpipe_biorep_label(label: str) -> str:
    """Sanitize Source Name or Sample Name for TMT BioReplicate (MSstatsTMT/FPAR). LFQ manifest uses numeric_biorep_for_subject instead."""
    if not label or not str(label).strip():
        return "1"
    return re.sub(r"[^A-Za-z0-9_]", "_", str(label).strip()) or "1"


def numeric_biorep_for_subject(subject_key: str, subject_to_id: dict, counter: list) -> str:
    """Assign stable numeric Bioreplicate IDs for FragPipe manifest (must be integer-like)."""
    key = (subject_key or "").strip()
    if not key:
        key = "__default__"
    if key not in subject_to_id:
        subject_to_id[key] = str(counter[0])
        counter[0] += 1
    return subject_to_id[key]
