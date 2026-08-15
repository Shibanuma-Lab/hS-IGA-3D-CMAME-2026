import os
import subprocess
import shutil
from const import simulation_params as sp
from utils.logger import logger
from utils.step2str import step2str


def run(step):
    """
    Execute SFEM solver for the given step
    
    Workflow:
    1. Copy inputfiles/step{XXXXX} to sfem_linear/example/
    2. Run solver in sfem_linear/example/step{XXXXX}/ with ../../bin/sfem_linear
    3. Move results back to circular_crack/results/step{XXXXX}/
    """
    logger.info(f"Running SFEM solver for step {step}")
    
    str_step = step2str(step)
    
    # Paths
    current_dir = os.getcwd()  # circular_crack/
    input_dir = os.path.join(current_dir, "inputfiles", f"step{str_step}")
    
    # SFEM paths (parent of circular_crack is project root)
    project_root = os.path.dirname(current_dir)
    sfem_dir = os.path.join(project_root, sp.REPO_NAME)
    sfem_example_dir = os.path.join(sfem_dir, "example")
    sfem_step_dir = os.path.join(sfem_example_dir, f"step{str_step}")
    sfem_executable = os.path.join(sfem_dir, "bin", "sfem_linear")
    
    # Results directory
    results_dir = os.path.join(current_dir, "results")
    results_step_dir = os.path.join(results_dir, f"step{str_step}")
    
    # Check if input directory exists
    if not os.path.exists(input_dir):
        logger.error(f"Input directory not found: {input_dir}")
        return False
    
    # Check if executable exists
    if not os.path.exists(sfem_executable):
        logger.error(f"SFEM executable not found: {sfem_executable}")
        logger.info("Please compile the SFEM solver first")
        return False
    
    # Create example directory if not exists
    os.makedirs(sfem_example_dir, exist_ok=True)
    
    # Step 1: Copy input directory to sfem_linear/example/
    logger.info(f"Copying {input_dir} to {sfem_example_dir}/")
    try:
        if os.path.exists(sfem_step_dir):
            shutil.rmtree(sfem_step_dir)
        shutil.copytree(input_dir, sfem_step_dir)
        logger.info(f"✓ Input files copied to {sfem_step_dir}")
    except Exception as e:
        logger.error(f"Failed to copy input files: {e}")
        return False
    
    # Step 2: Run solver in sfem_linear/example/step{XXXXX}/
    logger.info(f"Running solver in {sfem_step_dir}")
    original_dir = os.getcwd()
    
    try:
        os.chdir(sfem_step_dir)
        
        # Set OpenMP threads
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(sp.OPENMP)
        if 'SFEM_MASS_LUMPING_ALPHA' not in env and 'SFEM_MASS_LUMPING' not in env:
            env['SFEM_MASS_LUMPING_ALPHA'] = f"{sp.SFEM_MASS_LUMPING_ALPHA:.12g}"
        
        # Prepare command (relative path from step{XXXXX}/ to bin/sfem_linear)
        cmd = ['../../bin/sfem_linear', 'input.dat']
        
        logger.info(f"Executing: {' '.join(cmd)}")
        logger.info(f"Working directory: {os.getcwd()}")
        logger.info(f"OpenMP threads: {sp.OPENMP}")
        logger.info(
            "SFEM mass-lumping alpha: "
            f"{env.get('SFEM_MASS_LUMPING_ALPHA', 'legacy override')}"
        )
        
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
        else:
            # Run with terminal output
            result = subprocess.run(
                cmd,
                env=env,
                check=False
            )
        
        # Check return code
        if result.returncode == 0:
            logger.info(f"✓ SFEM solver completed successfully for step {step}")
            success = True
        else:
            logger.error(f"✗ SFEM solver failed with return code {result.returncode}")
            success = False
        
    except Exception as e:
        logger.error(f"Error executing SFEM solver: {e}")
        success = False
    
    finally:
        os.chdir(original_dir)
    
    # Step 3: Move results back to circular_crack/results/
    if success or os.path.exists(sfem_step_dir):
        logger.info(f"Moving results from {sfem_step_dir} to {results_step_dir}")
        try:
            # Create results directory if not exists
            os.makedirs(results_dir, exist_ok=True)
            
            # Remove old results if exists
            if os.path.exists(results_step_dir):
                shutil.rmtree(results_step_dir)
            
            # Move the entire step directory from example to results
            shutil.move(sfem_step_dir, results_step_dir)
            logger.info(f"✓ Results moved to {results_step_dir}")
            
        except Exception as e:
            logger.error(f"Failed to move results: {e}")
            # Try to at least copy if move fails
            try:
                if os.path.exists(results_step_dir):
                    shutil.rmtree(results_step_dir)
                shutil.copytree(sfem_step_dir, results_step_dir)
                shutil.rmtree(sfem_step_dir)
                logger.info(f"✓ Results copied to {results_step_dir}")
            except Exception as e2:
                logger.error(f"Failed to copy results: {e2}")
    
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
