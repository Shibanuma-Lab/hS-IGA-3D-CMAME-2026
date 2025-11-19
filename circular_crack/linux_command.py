import os
import subprocess
from const import simulation_params as sp
from utils.logger import logger
from utils.step2str import step2str


def run(step):
    """
    Execute SFEM solver for the given step
    """
    logger.info(f"Running SFEM solver for step {step}")
    
    str_step = step2str(step)
    input_dir = f"inputfiles/step{str_step}"
    
    # Check if input directory exists
    if not os.path.exists(input_dir):
        logger.error(f"Input directory not found: {input_dir}")
        return False
    
    # Path to SFEM executable
    sfem_executable = f"../../{sp.REPO_NAME}/bin/sfem_linear"
    
    # Check if executable exists
    if not os.path.exists(sfem_executable):
        logger.error(f"SFEM executable not found: {sfem_executable}")
        logger.info("Please compile the SFEM solver first")
        return False
    
    # Change to input directory
    original_dir = os.getcwd()
    os.chdir(input_dir)
    
    try:
        # Set OpenMP threads
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(sp.OPENMP)
        
        # Prepare command
        cmd = [sfem_executable, 'input.dat']
        
        logger.info(f"Executing: {' '.join(cmd)}")
        logger.info(f"OpenMP threads: {sp.OPENMP}")
        
        # Execute SFEM solver
        if sp.DOS_OPEN == 0:
            # Run without opening terminal
            with open('solver.log', 'w') as log_file:
                result = subprocess.run(
                    cmd,
                    env=env,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=False
                )
        elif sp.DOS_OPEN == 1:
            # Run and close terminal after completion
            result = subprocess.run(
                cmd,
                env=env,
                check=False
            )
        else:
            # Keep terminal open (sp.DOS_OPEN == 2)
            result = subprocess.run(
                cmd,
                env=env,
                check=False
            )
        
        # Check return code
        if result.returncode == 0:
            logger.info(f"SFEM solver completed successfully for step {step}")
            success = True
        else:
            logger.error(f"SFEM solver failed with return code {result.returncode}")
            success = False
        
    except Exception as e:
        logger.error(f"Error executing SFEM solver: {e}")
        success = False
    
    finally:
        # Return to original directory
        os.chdir(original_dir)
    
    return success


def run_parallel(step_list):
    """
    Run multiple steps in parallel (if applicable)
    Note: Typically serial for crack propagation
    """
    logger.info(f"Running SFEM solver for steps: {step_list}")
    
    results = []
    for step in step_list:
        success = run(step)
        results.append((step, success))
    
    return results


def check_solver_output(step):
    """
    Check if solver output files exist
    """
    str_step = step2str(step)
    input_dir = f"inputfiles/step{str_step}"
    
    required_files = ['delta_u.dat', 'stress.dat', 'strain.dat']
    
    missing_files = []
    for filename in required_files:
        filepath = os.path.join(input_dir, filename)
        if not os.path.exists(filepath):
            missing_files.append(filename)
    
    if missing_files:
        logger.warning(f"Missing output files for step {step}: {missing_files}")
        return False
    else:
        logger.info(f"All output files present for step {step}")
        return True


def extract_results(step):
    """
    Extract and organize results from solver output
    """
    logger.info(f"Extracting results for step {step}")
    
    str_step = step2str(step)
    input_dir = f"inputfiles/step{str_step}"
    output_dir = f"{sp.OUTPUT_FOLDER}/step{str_step}"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Copy important result files
    import shutil
    
    files_to_copy = [
        'delta_u.dat',
        'velocity.dat', 
        'acceleration.dat',
        'stress.dat',
        'strain.dat',
        'reaction.dat'
    ]
    
    for filename in files_to_copy:
        src = os.path.join(input_dir, filename)
        dst = os.path.join(output_dir, filename)
        
        if os.path.exists(src):
            shutil.copy2(src, dst)
            logger.info(f"Copied {filename} to output folder")
        else:
            logger.warning(f"File not found: {src}")
    
    logger.info(f"Results extracted to {output_dir}")


def cleanup_temp_files(step):
    """
    Clean up temporary files after successful run
    """
    logger.info(f"Cleaning up temporary files for step {step}")
    
    str_step = step2str(step)
    input_dir = f"inputfiles/step{str_step}"
    
    # List of temporary files to remove
    temp_files = [
        'solver.tmp',
        'matrix.tmp',
        'debug.log'
    ]
    
    for filename in temp_files:
        filepath = os.path.join(input_dir, filename)
        if os.path.exists(filepath):
            os.remove(filepath)
            logger.info(f"Removed temporary file: {filename}")
