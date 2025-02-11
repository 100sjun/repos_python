# Modules
import subprocess

# compile zdplaskin
def run_prep(inp_path):
    try:
        process = subprocess.Popen(
            'preprocessor.exe',
            stdout = subprocess.DEVNULL,        # ignore outputs
            stderr = subprocess.DEVNULL,        # ignore errors
            stdin = subprocess.PIPE,            # recognize input
            universal_newlines=True
        )
        
        process.stdin.write(inp_path)
        process.stdin.flush()                   # send a data

        while process.poll() is None:           # check the program state, if None, program is still in the run
            process.stdin.write('.\n')
            process.stdin.flush()
    except:
        pass
    print('check the run of preprocessor')
    return process

# Compile exe
def compile_zdp(name):
    compile_command = [
        'gfortran', '-o', name, 'dvode_f90_m.F90', 'zdplaskin_m.F90',
        'run.F90', 'bolsig_x86_64_g.dll'
    ]
    
    try:
        subprocess.run(compile_command)
        
    except:
        pass
    print('check the compiler')

# Run exe
def run_exe(exe_path):
    try:
        process = subprocess.Popen(
            exe_path,
            stdout = subprocess.PIPE, # Read standard outputs
            stderr = subprocess.PIPE, # Read standard errors
            stdin = subprocess.PIPE,            # recognize input
            universal_newlines=True,  # outputs to str variables
            bufsize = 1               # control the size of buffer
        )

        log_flag = False             # The flag for starting log after "Caculation Start!!"
        while True:
            output = process.stdout.readline()
            if not output:
                break
            if "Calculation Start" in output:
                log_flag = True

            if log_flag:
                print(f'\r{output.strip()}           ',end='',flush=True)
            if "Press ENTER to continue" in output:
                process.stdin.write('\n')
                process.stdin.flush()  # Flush the buffer to ensure the input is sent
                
            if "PRESS ENTER TO EXIT" in output:
                process.kill()        # forced shutdown
                break
            if "WARNING: BOLSIG+ convergence failed" in output:
                process.stdin.write('\n')
                process.stdin.flush()  # Flush the buffer to ensure the input is sent
    except:
        pass
    return process

inp_path = 'kinet.inp'
exe_path = 'run.exe'
prep_process = run_prep(inp_path)
prep_process.wait()
compile_zdp(exe_path)
run_exe(exe_path)