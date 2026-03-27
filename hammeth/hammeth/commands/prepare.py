import subprocess
import sys
from pathlib import Path


def run_prepare(args):

    # 定位 scripts/prepareData.py
    script_path = Path(__file__).resolve().parents[1] / "scripts" / "prepareData.py"

    cmd = [
        sys.executable,
        str(script_path),
        "-i", args.i,
        "-m", str(args.m),
        "-r", args.r,
        "-cp", args.cp,
        "-g", args.g,
        "-cl",
    ] + args.cl

    if args.o:
        cmd += ["-o", args.o]

    print("Running command:")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)
