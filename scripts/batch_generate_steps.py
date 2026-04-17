#!/usr/bin/env python3
"""Batch-generate STEP files for all MAS examples using MVB++ with MKF enrichment."""

import sys
import os
import subprocess
import json
from pathlib import Path

# Paths
MVBPP_BUILD = Path("/home/alf/OpenMagnetics/MVB++/build")
MVBPP_GENERATOR = MVBPP_BUILD / "mvbpp_step_generator"
MAS_EXAMPLES = Path("/home/alf/OpenMagnetics/MAS/examples")
DEFAULT_OUTPUT = Path("/home/alf/OpenMagnetics/MVB++/mas_example_steps")


def run_mvbpp(input_file, output_file, use_mkf=True):
    """Run MVB++ step generator on a single file."""
    cmd = [
        str(MVBPP_GENERATOR),
        str(input_file),
        "-o", str(output_file)
    ]
    if not use_mkf:
        cmd.append("--no-mkf")
    
    env = {**os.environ, "LD_LIBRARY_PATH": str(MVBPP_BUILD / "occt-install/lib")}
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, env=env)
        return result.returncode == 0, result.stderr
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)


def main():
    parser = argparse.ArgumentParser(description="Batch-generate STEP files for MAS examples")
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT, help="Output directory")
    parser.add_argument("--no-mkf", action="store_true", help="Skip MKF enrichment")
    parser.add_argument("examples", nargs="*", help="Specific example files to process")
    args = parser.parse_args()
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {args.output}")
    
    # Get list of examples
    if args.examples:
        examples = [Path(e) for e in args.examples]
    else:
        examples = sorted(MAS_EXAMPLES.glob("*.json"))
    
    print(f"\nProcessing {len(examples)} examples with MKF mas_autocomplete...\n")
    
    success_count = 0
    failed = []
    
    for i, ex in enumerate(examples, 1):
        print(f"[{i}/{len(examples)}] {ex.stem}...", end=" ", flush=True)
        
        out_file = args.output / f"{ex.stem}_mvbpp.step"
        
        # Try with MKF (mas_autocomplete)
        if not args.no_mkf:
            success, err = run_mvbpp(ex, out_file, use_mkf=True)
            if success:
                print("✓ (MKF)")
                success_count += 1
                continue
            else:
                # If MKF fails, try without
                print(f"MKF failed, trying without...", end=" ")
        
        # Try without MKF
        success, err = run_mvbpp(ex, out_file, use_mkf=False)
        if success:
            print("✓ (raw)")
            success_count += 1
        else:
            print(f"✗ Failed: {err[:60]}")
            failed.append((ex.stem, err))
    
    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY: {success_count}/{len(examples)} succeeded")
    if failed:
        print(f"\nFailed examples ({len(failed)}):")
        for name, err in failed[:10]:  # Show first 10
            print(f"  - {name}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more")
    print(f"\nSTEP files saved to: {args.output}")
    
    return 0 if len(failed) == 0 else 1


if __name__ == "__main__":
    import argparse
    sys.exit(main())
