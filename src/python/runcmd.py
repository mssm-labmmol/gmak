import subprocess

def run(cmd):
    # Identify pipes with 'echo'
    piped = cmd.split('|')
    if len(piped) == 1:
        # No pipes!
        return subprocess.run(cmd.split()).check_returncode()
    else:
        # Execute command before pipe, storing stdout.
        command2input = subprocess.check_output(piped[0].split())
        # Open process.
        p = subprocess.Popen(piped[1].split(), stdin=subprocess.PIPE)
        # Communicate stdin.
        p.communicate(input=command2input)
