import subprocess

class ConfigVariables:

    # get gmx path from 'which gmx'
    gmx = subprocess.check_output(['which', 'gmx']).strip().decode('utf-8')
    # default temporary directory to create files and directories inside
    tmpdir = '/tmp'
    

    
