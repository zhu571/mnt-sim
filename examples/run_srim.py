#!/usr/bin/env python
"""
Run SRIM to calculate stopping power for 243U in 238U target.
"""
import subprocess, os, shutil, glob, time

srim_dir = "D:/software/SRIM"
sr_module_dir = os.path.join(srim_dir, "SR Module")

# Create SR.IN for 243U in Uranium
sr_in_content = """---Stopping/Range Input Data (Number-format: Period = Decimal Point)
---Output File Name
"243U in Uranium - Full Range"
---Ion(Z), Ion Mass(u)
92   243.0
---Target Data: (Solid=0,Gas=1), Density(g/cm3), Compound Corr.
0    18.95   1.0
---Number of Target Elements
 1 
---Target Elements: (Z), Target name, Stoich, Target Mass(u)
92   "Uranium"               1             238.0
---Output Stopping Units (1-8)
 6
---Ion Energy : E-Min(keV), E-Max(keV)
 10    2000000
"""

# Write SR.IN to SR Module directory
with open(os.path.join(sr_module_dir, "SR.IN"), "w") as f:
    f.write(sr_in_content)

# Also write to main SRIM dir
with open(os.path.join(srim_dir, "SR.IN"), "w") as f:
    f.write(sr_in_content)

print("SR.IN written. Running SRModule.exe...")

# Run SRModule.exe
result = subprocess.run(
    [os.path.join(sr_module_dir, "SRModule.exe")],
    cwd=sr_module_dir,
    capture_output=True, text=True, timeout=120
)
print(f"Exit code: {result.returncode}")
if result.stdout: print("STDOUT:", result.stdout[:500])
if result.stderr: print("STDERR:", result.stderr[:500])

# Wait a moment for file to be written
time.sleep(1)

# Find output files
output_dir = os.path.join(srim_dir, "SR Outputs")
if os.path.exists(output_dir):
    files = [f for f in os.listdir(output_dir) if '243U' in f or 'Uranium' in f]
    print(f"\nOutput files found: {files}")
    for fname in files:
        fpath = os.path.join(output_dir, fname)
        with open(fpath) as f:
            content = f.read(2000)
        print(f"\n--- {fname} ---")
        print(content)

# Also search for any new files
print("\nSearching for output files...")
for root, dirs, files in os.walk(srim_dir):
    for f in files:
        if '243U' in f or 'Uranium' in f or f.endswith('.txt'):
            fpath = os.path.join(root, f)
            mtime = os.path.getmtime(fpath)
            if time.time() - mtime < 300:  # created in last 5 min
                print(f"  {fpath} ({os.path.getsize(fpath)} bytes)")
