import subprocess
import os
from gmak.config import ConfigVariables
import re

no_backup_environ = dict(os.environ, GMX_MAXBACKUP="-1")

def run(cmd):
    # Identify pipes with 'echo'
    piped = cmd.split('|')
    if len(piped) == 1:
        # No pipes!
        return subprocess.run(cmd.split(), env=no_backup_environ).check_returncode()
    else:
        # Execute command before pipe, storing stdout.
        command2input = subprocess.check_output(piped[0].split())
        # Open process.
        p = subprocess.Popen(piped[1].split(), stdin=subprocess.PIPE,
                             env=no_backup_environ)
        # Communicate stdin.
        p.communicate(input=command2input)


def gmx_insert_molecules(args_list):
    default_try = '1000'
    from gmak.configurations import ConfigurationError
    # number of molecules requested
    nmol = int(args_list[args_list.index("-nmol") + 1])
    # set try if not given in args_list
    if '-try' not in args_list:
        args_list += ['-try', default_try]
    # save stdout to check if requested number of molecules was
    # reached
    log = subprocess.check_output(
        [ConfigVariables.gmx, "insert-molecules"] + args_list,
        env=no_backup_environ,
        stderr=subprocess.STDOUT).decode('utf-8')
    # check number of added molecules
    m = re.search('Added ([0-9])+ molecules', log)
    nmol_added = int(m.group(1))
    if nmol != nmol_added:
        raise ConfigurationError(f"Requested {nmol} molecules but "
                                 "gmx insert-molecules was able "
                                 f"to insert only {nmol_added}.")
    
    
