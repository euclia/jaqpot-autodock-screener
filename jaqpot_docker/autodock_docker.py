import jaqpotpy as jp

# import deepchem as jp
# from deepchem.dock.pose_generation import VinaPoseGenerator
# from deepchem.dock import Docker

from jaqpotpy.docking.pose_generation import VinaPoseGenerator
from jaqpotpy.docking import Docker


class SimpleDocker:

    def __init__(self, out_dir, exhaustiveness, num_modes, calc_charges:bool = False, add_hydrogens:bool = False):
        self.out_dir = out_dir
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.calc_charges = calc_charges
        self.add_hydrogens = add_hydrogens

    def dock_ligand(self, protein_file, ligand_file):
        vpg = VinaPoseGenerator(calc_charges=self.calc_charges, add_hydrogens=self.add_hydrogens)
        docker = Docker(vpg)
        docked_outputs = docker.dock((protein_file, ligand_file),
                                         exhaustiveness=self.exhaustiveness,
                                         num_modes=self.num_modes,
                                         out_dir=self.out_dir,
                                         use_pose_generator_scores=True)
        docked_outputs = list(docked_outputs)
        return docked_outputs

    def dock_async(self, protein_file, ligand_file):
        vpg = VinaPoseGenerator(calc_charges=self.calc_charges, add_hydrogens=self.add_hydrogens)
        docker = Docker(vpg)
        docked_outputs = docker.dock((protein_file, ligand_file),
                                         exhaustiveness=self.exhaustiveness,
                                         num_modes=self.num_modes,
                                         out_dir=self.out_dir,
                                         use_pose_generator_scores=True)
        return docked_outputs

    def dock_ligand_in_pocket(self, protein_file, ligand_file, centroid, box_dims):
        vpg = VinaPoseGenerator(calc_charges=self.calc_charges, add_hydrogens=self.add_hydrogens)
        docker = Docker(vpg)
        docked_outputs = docker.dock((protein_file, ligand_file),
                                     centroid=centroid,
                                     box_dims=box_dims,
                                     exhaustiveness=self.exhaustiveness,
                                     num_modes=self.num_modes,
                                     out_dir=self.out_dir,
                                     use_pose_generator_scores=True)
        return docked_outputs

    def dock_ligand_in_pocket_sync(self, protein_file, ligand_file, centroid, box_dims):
        vpg = VinaPoseGenerator(calc_charges=self.calc_charges, add_hydrogens=self.add_hydrogens)
        docker = Docker(vpg)
        docked_outputs = docker.dock((protein_file, ligand_file),
                                     centroid=centroid,
                                     box_dims=box_dims,
                                     exhaustiveness=self.exhaustiveness,
                                     num_modes=self.num_modes,
                                     out_dir=self.out_dir,
                                     use_pose_generator_scores=True)
        docked_outputs = list(docked_outputs)
        return docked_outputs


    def extract_site(self, protein_file, ligand_file):
        active_site_box, active_site_coords = dock.binding_pocket.extract_active_site(protein_file, ligand_file)
        return active_site_box, active_site_coords