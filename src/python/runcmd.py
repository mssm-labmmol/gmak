import subprocess

def run(cmd):
    return subprocess.run(cmd.split()).check_returncode()
