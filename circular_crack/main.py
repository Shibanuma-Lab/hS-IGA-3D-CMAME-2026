import os
import argparse
import numpy as np
from const import const_local_mesh, const_global_mesh, simulation_params as sp
from const import const_jintegral as cji
from initial import initial
import global_mesh
import local_mesh
import boundary
import input_generator
from scripts import linux_command
from jintegral import JIntegralCalculator
from utils.logger import logger
from utils.step2str import step2str

def makemeshs(step, REstart):
    logger.info(os.getcwd())
    str_step = step2str(step)
    
    # Create necessary directories
    os.makedirs("inputfiles", exist_ok=True)
    os.makedirs(f"inputfiles/step{str_step}", exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    
    os.chdir(f"inputfiles/step{str_step}")
    logger.info(f"step: {step} :: Generate Global Mesh")
    g = global_mesh.GlobalMesh(step)
    g.make_global_mesh()
    g.generate()

    logger.info(f"step: {step} :: Generate Local Mesh")
    l = local_mesh.LocalMesh(step)
    l.make_local_mesh()
    l.generate()

    logger.info(f"step: {step} :: Generate Boundary")
    b = boundary.Boundary(l, g)
    logger.info(f"step: {step} :: Define Boundary and Load")
    b.define_boundary(l, g)
    b.generate()

    logger.info(f"step: {step} :: Generate Input File")
    input_generator.generate(step, REstart)
    
    # Generate virtual mesh if needed
    input_generator.generate_virtual_mesh()
    
    if REstart == 1:
        logger.info(f"step: {step} :: Initialize from Previous Step")
        init = initial(step, l, g)
    
    os.chdir("../../")
    return l, g


# makemeshs(200)
# linux_command.run(200)

def main():
    """
    Main function with command line argument support
    
    Usage examples:
        # Run default steps from config
        python3 main.py
        
        # Generate mesh only for step 2
        python3 main.py --step_start 2 --step_end 3 --meshonly
        
        # Run single step with restart
        python3 main.py --particular --step_start 5
        
        # Run steps 0-10 with solver
        python3 main.py --step_start 0 --step_end 10
        
        # Calculate J-integral and K factors from results
        python3 main.py --is_K --step_start 1 --step_end 101
        
        # Calculate J-integral with custom parameters
        python3 main.py --is_K --step_start 50 --step_end 101 --Rj0 1.5 --Rj1 2.0
        
        # Debug mode (more logging)
        python3 main.py --debugmode --step_start 2 --step_end 3
    """
    argparser = argparse.ArgumentParser(
        description="S-IGA Circular Crack Simulation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run simulation
  python3 main.py --step_start 2 --step_end 3 --meshonly
  python3 main.py --particular --step_start 5
  python3 main.py --step_start 0 --step_end 10
  
  # Post-process results (calculate J-integral and K factors)
  python3 main.py --is_K --step_start 1 --step_end 101
  python3 main.py --is_K --step_start 50 --step_end 101 --Rj0 1.5 --Rj1 2.0
        """
    )
    
    # Step range control
    argparser.add_argument("--step_start", type=int, default=sp.step_start,
                          help=f"Starting step number (default: {sp.step_start})")
    argparser.add_argument("--step_end", type=int, default=sp.step_end,
                          help=f"Ending step number (default: {sp.step_end})")
    
    # Execution modes
    argparser.add_argument("--meshonly", action="store_true", default=False,
                          help="Only generate mesh and input files, skip solver")
    argparser.add_argument("--solveronly", action="store_true", default=False,
                          help="Only run solver (assumes input files exist)")
    argparser.add_argument("--particular", action="store_true", default=False,
                          help="Run single particular step (only step_start)")
    
    # Restart control
    argparser.add_argument("--restart", type=int, default=None,
                          help="Force REstart value (0: fresh start, 1: restart from previous)")
    argparser.add_argument("--no_restart", action="store_true", default=False,
                          help="Force fresh start for all steps (REstart=0)")
    
    # Debug and cleanup
    argparser.add_argument("--debugmode", action="store_true", default=False,
                          help="Enable debug mode (more verbose logging)")
    argparser.add_argument("--delete", action="store_true", default=False,
                          help="Delete all files in inputfiles/ before running")
    
    # J-integral and K calculation
    argparser.add_argument("--is_K", action="store_true", default=False,
                          help="Calculate J-integral and stress intensity factors from results")
    argparser.add_argument("--Rj0", type=float, default=None,
                          help=f"Inner radius for J-integral domain (default: {cji.Rj0})")
    argparser.add_argument("--Rj1", type=float, default=None,
                          help=f"Outer radius for J-integral domain (default: {cji.Rj1})")
    argparser.add_argument("--Wj0", type=float, default=None,
                          help=f"Inner width for J-integral domain (default: {cji.Wj0})")
    argparser.add_argument("--Wj1", type=float, default=None,
                          help=f"Outer width for J-integral domain (default: {cji.Wj1})")
    argparser.add_argument("--velocity", type=float, default=None,
                          help=f"Crack velocity for J-integral [m/s] (default: {sp.V})")
    argparser.add_argument("--output_K", type=str, default=None,
                          help=f"Output file for J-integral results (default: {cji.output_file})")
    
    args = argparser.parse_args()
    
    # Log arguments with actual default values displayed
    logger.info("="*60)
    logger.info("S-IGA Circular Crack Simulation")
    logger.info("="*60)
    logger.info(f"Command line arguments:")
    
    # Display with actual defaults instead of None for better clarity
    display_args = vars(args).copy()
    if args.is_K:
        # Show actual defaults for J-integral parameters when in K calculation mode
        if display_args['Rj0'] is None:
            display_args['Rj0'] = f"{cji.Rj0} (default)"
        if display_args['Rj1'] is None:
            display_args['Rj1'] = f"{cji.Rj1} (default)"
        if display_args['Wj0'] is None:
            display_args['Wj0'] = f"{cji.Wj0} (default)"
        if display_args['Wj1'] is None:
            display_args['Wj1'] = f"{cji.Wj1} (default)"
        if display_args['velocity'] is None:
            display_args['velocity'] = f"{sp.V} (default)"
        if display_args['output_K'] is None:
            display_args['output_K'] = f"{cji.output_file} (default)"
    
    for arg, value in display_args.items():
        logger.info(f"  {arg}: {value}")
    logger.info("="*60)
    
    # Handle delete flag
    if args.delete:
        import shutil
        logger.warning("DELETE mode enabled - removing inputfiles/ directory")
        if os.path.exists("inputfiles"):
            shutil.rmtree("inputfiles")
            logger.info("inputfiles/ directory deleted")
        else:
            logger.info("inputfiles/ directory does not exist")
    
    # Handle particular mode
    if args.particular:
        args.step_end = args.step_start + 1
        logger.info(f"PARTICULAR mode: running only step {args.step_start}")
    
    # Handle J-integral calculation mode
    if args.is_K:
        logger.info("="*60)
        logger.info("J-INTEGRAL AND STRESS INTENSITY FACTOR CALCULATION MODE")
        logger.info("="*60)
        
        # Use default values from const_jintegral if not specified
        step_start = args.step_start if args.step_start != sp.step_start else cji.step_start
        step_end = args.step_end if args.step_end != sp.step_end else cji.step_end
        Rj0 = args.Rj0 if args.Rj0 is not None else cji.Rj0
        Rj1 = args.Rj1 if args.Rj1 is not None else cji.Rj1
        Wj0 = args.Wj0 if args.Wj0 is not None else cji.Wj0
        Wj1 = args.Wj1 if args.Wj1 is not None else cji.Wj1
        velocity = args.velocity if args.velocity is not None else sp.V
        output_K = args.output_K if args.output_K is not None else cji.output_file
        
        logger.info(f"Reading results from steps {step_start} to {step_end}")
        logger.info(f"J-domain parameters: Rj=[{Rj0}, {Rj1}], Wj=[{Wj0}, {Wj1}]")
        logger.info(f"Crack velocity: {velocity} m/s")
        logger.info(f"Output file: {output_K}")
        
        try:
            calculator = JIntegralCalculator(
                step_start=step_start,
                step_end=step_end,
                Rj0=Rj0,
                Rj1=Rj1,
                Wj0=Wj0,
                Wj1=Wj1,
                v=velocity
            )
            calculator.run(output_file=output_K)
        except Exception as e:
            logger.error(f"J-integral calculation failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
        
        return
    
    # Determine step range
    step_start = args.step_start
    step_end = args.step_end
    
    if step_start >= step_end:
        logger.error(f"Invalid step range: start={step_start}, end={step_end}")
        return
    
    logger.info(f"Running steps {step_start} to {step_end-1}")
    
    # Main execution loop
    step_list = range(step_start, step_end)
    for step in step_list:
        logger.info("")
        logger.info("="*60)
        logger.info(f"STEP {step}")
        logger.info("="*60)
        
        # Determine REstart value
        if args.restart is not None:
            # User explicitly specified REstart value
            REstart = args.restart
        elif args.no_restart:
            # User wants to force fresh start
            REstart = 0
        else:
            # Default behavior:
            # - If step_start != 0, we're resuming from middle, so REstart=1 for first step
            # - For subsequent steps, always REstart=1
            if step == step_start:
                # First step in this run
                if step_start == 0:
                    REstart = 0  # True fresh start from step 0
                else:
                    REstart = 1  # Resuming from middle, need initial conditions
            else:
                REstart = 1  # All subsequent steps restart from previous
        
        logger.info(f"REstart: {REstart} (step={step}, step_start={step_start})")
        
        # Generate mesh and input files (unless solveronly mode)
        if not args.solveronly:
            logger.info(f"Generating mesh and input files for step {step}")
            l, g = makemeshs(step, REstart)
        else:
            logger.info(f"Skipping mesh generation (solveronly mode)")
        
        # Run solver (unless meshonly mode)
        if not args.meshonly:
            logger.info(f"Running solver for step {step}")
            linux_command.run(step)
        else:
            logger.info(f"Skipping solver (meshonly mode)")
    
    logger.info("")
    logger.info("="*60)
    logger.info("ALL STEPS COMPLETED")
    logger.info("="*60)


if __name__ == "__main__":
    main()