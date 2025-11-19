import os
import numpy as np
from const import const_local_mesh, const_global_mesh, simulation_params as sp
from initial import initial
import global_mesh
import local_mesh
import boundary
import input_generator
import linux_command
# from jintegral import jintegral
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
    step_start = sp.step_start
    step_end = sp.step_end

    step_list = range(step_start, step_end)
    for step in step_list:
        REstart = 0 if step == step_start else 1  # 1:restart from dynamic analysis
        logger.info(f"Step {step}, REstart: {REstart}")
        l, g = makemeshs(step, REstart)
        linux_command.run(step)


if __name__ == "__main__":
    main()