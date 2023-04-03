import argparse
from minio import Minio
from minio.error import S3Error
import gzip
from rdkit import Chem
import pandas as pd
import ast
from jaqpot_docker import SimpleDocker


def create_minio_client(minio_service, minio_username, minio_password, secure):
    print("Creating s3 client")
    global minio_client
    global bucket_name
    if secure == "False":
        minio_client = Minio(
                minio_service,
                access_key=minio_username,
                secret_key=minio_password,
            secure=False
            )
    if secure == "True":
        minio_client = Minio(
                minio_service,
                access_key=minio_username,
                secret_key=minio_password,
            secure=True
            )
    bucket_name = "zinc012023"
    return minio_client, bucket_name


def download_pdb(pdb_file):
    print("Downloading pdb")
    bucket = pdb_file.split("/")[0]
    pdb_file_name = pdb_file.split("/")[1]
    global protein_pdb
    protein_pdb = "./" + pdb_file_name
    minio_client.fget_object(bucket, pdb_file_name, protein_pdb)


def download_sdf(molecule_sdf):
    print("Downloading sdf")
    bucket = molecule_sdf.split("/")[0]
    molecule_sdf_file = molecule_sdf.split("/")[1]
    global molecule_sdf_local
    molecule_sdf_local = "./" + molecule_sdf_file
    minio_client.fget_object(bucket, molecule_sdf_file, molecule_sdf_local)
    f = gzip.open(molecule_sdf_local,'rb')
    file_content=f.read()
    global sdf_file
    sdf_file = molecule_sdf_local[:len(molecule_sdf_local)-3]
    with open(sdf_file, "wb") as sdf:
        sdf.write(file_content)


def create_simple_docker(exhaustiv, num_modes, calc_charges=False, add_hydrogenes=False,  output="./tmp"):
    print("Creating docker")
    global docker
    if calc_charges == "True":
        calc_charges = True
    else:
        calc_charges = False
    if add_hydrogenes == "True":
        add_hydrogenes = True
    else:
        add_hydrogenes = False
    docker = SimpleDocker(exhaustiveness = int(exhaustiv)
                        , num_modes = int(num_modes)
                        , add_hydrogens = add_hydrogenes
                        , calc_charges = calc_charges
                        , out_dir = output)
    return docker


def start_docking(task_bucket, centroid, box_dims):
    print("Starting docking")
    centroid  = ast.literal_eval(centroid)
    box_dims = ast.literal_eval(box_dims)
    df = pd.DataFrame()
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        run = True
        if mol is not None:
            uniqueId = mol.GetProp('_Name')
            ligand_file = "./tmp/" + uniqueId + ".sdf"
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in ['Cl', "Br"]:
                    run = False
                    continue
            if run is True:
                with Chem.SDWriter(ligand_file) as w:
                    w.write(mol)
                docked = docker.dock_ligand_in_pocket_sync( protein_pdb
                                                            , ligand_file
                                                            , centroid=(centroid[0], centroid[1], centroid[2])
                                                            , box_dims=(box_dims[0], box_dims[1], box_dims[2]))
                df = df.append({"Docking_molecule": Chem.MolToSmiles(mol, kekuleSmiles=True)
                                        , "Unique_ID": str(uniqueId)
                                        , "Docked_value_bfe": docked[0][1]
                                        , "Ligand_File_Docked": uniqueId + "_docked.pdbqt"
                                        }
                                    , ignore_index=True
                                    )
                df.to_csv('./results/results.csv', index=False)

def main(args):
    task_id = args.task_id
    molecule_sdf_file = args.molecule_sdf_file
    molecule_smi_file = args.molecule_smi_file
    protein_pdb = args.protein_pdb
    centroid = args.centroid
    box_dims = args.box_dims
    exhaustiveness = args.exhaustiveness
    add_hydrogenes= args.add_hydrogens
    calc_charges = args.calc_charges
    num_modes = args.num_modes
    s3_service = args.s3_service
    s3_secure = args.s3_secure
    s3_username = args.s3_username
    s3_password = args.s3_password
    create_minio_client(s3_service, s3_username, s3_password, s3_secure)
    download_pdb(protein_pdb)
    download_sdf(molecule_sdf_file)
    create_simple_docker(exhaustiv=exhaustiveness, num_modes=num_modes, calc_charges=calc_charges, add_hydrogenes=add_hydrogenes)
    start_docking(task_bucket=task_id, centroid=centroid, box_dims=box_dims)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Jaqpot Autodock Vina screener",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t", "--task-id", dest="task_id", help="Task id. The task also is the bucket on the s3")
    parser.add_argument("-m", "--molecule-sdf", dest="molecule_sdf_file", help="The molecular sdf file and the bucket")
    parser.add_argument("-ms", "--molecule-smi", dest="molecule_smi_file", help="The molecular smi ( smiles ) file and the bucket")
    parser.add_argument("-pdb", "--pdb", dest="protein_pdb", help="The protein pdb file to dock and the bucket")
    parser.add_argument("-centroid", dest="centroid", help="The protein box centroid")
    parser.add_argument("-box_dims", dest="box_dims",  help="The protein box dimensions")
    parser.add_argument("-calc_charges", dest="calc_charges", default=False, help="Add charges on pdb")
    parser.add_argument("-add_hydrogens", dest="add_hydrogens", default=False, help="Add hydrogenes on pdb")
    parser.add_argument("-exhaust", dest="exhaustiveness", help="The exhaustiveness of the Autodock vina screener")
    parser.add_argument("-poses", dest="num_modes", help="The number of the poses to be generated")
    parser.add_argument("-s3", dest="s3_service", default="minio-service.dbs:9000", help="The s3 service")
    parser.add_argument("-s3s", dest="s3_secure", default=False, help="http/https s3 service")
    parser.add_argument("-usn", dest="s3_username", default="jaqpot", help="The s3 service username")
    parser.add_argument("-pass", dest="s3_password", default="jaqpot@minio", help="The s3 service password")
    args = parser.parse_args()
    config = vars(args)
    main(args)



#python3 main.py -centroid "[0.0, 0.0, 0.0]" -t "sample_dock" -m "zinc012023/3D_BA_AAML_BAAAML.xaa.sdf.gz" -pdb "sampledocking/AF-Q7PRQ0-F1-model_v4.pdb" -box_dims "[20, 20, 20]"  -calc_charges True -add_hydrogens True -s3 "localhost:9000" -s3s "False" -usn "minioadmin" -pass "minioadmin" -exhaust "1" -poses "1"