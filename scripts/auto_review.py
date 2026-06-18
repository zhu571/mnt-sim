#!/usr/bin/env python
"""
Scan mnt-sim for recent file changes and create Kanban review cards.
Designed to run as a cron job for adversarial review.
"""
import os
import sys
import time
import subprocess
import json

PROJECT = "D:/work/agent work/mnt-sim"
STATE_FILE = os.path.join(os.path.dirname(__file__), ".review_state.json")

# Track which files to review
WATCHED = [
    "mnt_sim/cross_section/imqmd_like.py",
    "mnt_sim/cross_section/imqmd_mc.py",
    "mnt_sim/cross_section/dns.py",
    "mnt_sim/cross_section/grazing.py",
    "mnt_sim/transport/gas_cell.py",
    "mnt_sim/transport/monte_carlo.py",
    "mnt_sim/transport/stopping.py",
]


def load_state():
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE) as f:
            return json.load(f)
    return {"last_check": 0, "reviewed_versions": {}}


def save_state(state):
    with open(STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)


def get_mtime(path):
    try:
        return os.path.getmtime(os.path.join(PROJECT, path))
    except OSError:
        return 0


def changed_files(state):
    changed = []
    for rel in WATCHED:
        mtime = get_mtime(rel)
        prev = state["reviewed_versions"].get(rel, 0)
        if mtime > prev:
            changed.append((rel, mtime))
    return changed


def create_review_card(file_list):
    files_str = "\n".join(f"  - D:/work/agent work/mnt-sim/{f}" for f, _ in file_list)
    body = f"""## 自动触发：文件变更审查

以下文件最近被修改，需要 adversarial review：

{files_str}

请重点审查：
1. 运动学公式是否正确
2. SRIM 数据加载和处理是否正确
3. 代码健壮性（路径、错误处理）
4. 可复现性（随机种子）
"""

    title = f"审查 MNT-SIM: {', '.join(f.split('/')[-1] for f, _ in file_list[:3])}"
    if len(file_list) > 3:
        title += f" +{len(file_list)-3} more"

    cmd = [
        "hermes", "kanban", "create",
        title,
        "--assignee", "auditor",
        "--body", body,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    return result.stdout, result.stderr


def main():
    state = load_state()
    state["last_check"] = time.time()

    changed = changed_files(state)
    if not changed:
        print("No changes detected.")
        save_state(state)
        return

    print(f"Changed files: {[f for f, _ in changed]}")
    out, err = create_review_card(changed)
    if err:
        print(f"Error: {err}")
    else:
        print(f"Review card created: {out.strip()}")

    # Update review state
    for rel, mtime in changed:
        state["reviewed_versions"][rel] = mtime
    save_state(state)


if __name__ == "__main__":
    main()
