"""Resolve logical frozen inputs to the main/supplement/shared input folders."""
from pathlib import Path
import json
def resolve(root,logical):
    root=Path(root);mapping=json.loads((root/'metadata/path_map.json').read_text(encoding='utf8'))
    return root/mapping[logical]
