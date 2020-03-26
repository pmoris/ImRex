""" Preprocess ppi files. """
from itertools import chain

import pandas as pd
from requests_futures.sessions import (
    FuturesSession,
    ThreadPoolExecutor,
)
from tqdm import tqdm

import src.bacli as bacli

PPI_NAME1 = "Unique identifier for interactor A"
PPI_NAME2 = "Unique identifier for interactor B"
PPI_NAMES = [PPI_NAME1, PPI_NAME2]

SEQUENCES_PATH = "PPI_sequences.csv"


@bacli.command  # noqa: C901
def ppi(
    inpath1: str,
    outpath1: str,
    inpath2: str = None,
    outpath2: str = None,
    sequences_path: str = SEQUENCES_PATH,
):
    def get_accession(identifier):
        if not identifier.startswith("uniprotkb:"):
            return None
        accession = identifier.split(":")[1]
        # DO NOT DROP PART AFTER dash (-)
        # accession = accession.split('-')[0]
        return accession

    def get_url(identifier):
        accession = get_accession(identifier)
        url = f"https://www.ebi.ac.uk/proteins/api/proteins/{accession}"
        return url

    class PpiDataset(object):
        def __init__(self, inpath, outpath):
            self.inpath = inpath
            self.outpath = outpath
            self.data = ppi_trim(pd.read_csv(self.inpath, sep="\t"))

        def get_identifiers(self):
            return set(chain(self.data[PPI_NAME1], self.data[PPI_NAME2]))

    datasets = [PpiDataset(inpath1, outpath1)]
    if inpath2 is not None:
        datasets.append(PpiDataset(inpath2, outpath2))

    identifiers = set(
        identifier for dataset in datasets for identifier in dataset.get_identifiers()
    )
    print("Unique identifiers: ", len(identifiers))

    sequences = list()
    failed = list()

    def get_response_handler(identifier):
        def handler(response, **kwargs):
            if not response.ok:
                tqdm.write(f"Failed request: {identifier}")
                tqdm.write("\t" + response.text)
                tqdm.write("removing it later")
                failed.append(identifier)
                # for dataset in datasets:
                #     for name in PPI_NAMES:
                #     dataset.data = dataset.data.drop(dataset.data[(dataset.data[PPI_NAME1] == identifier) | (dataset.data[PPI_NAME2] == identifier)].index)
                # dataset.data = dataset.data[dataset.data[name].str.equals(identifier)]
                return
            json_data = response.json()
            sequence = json_data["sequence"]["sequence"]

            sequences.append((identifier, sequence))

        return handler

    session = FuturesSession(executor=ThreadPoolExecutor(max_workers=20))

    def get_request(identifier):
        url = get_url(identifier)
        return session.get(url, hooks={"response": get_response_handler(identifier)})

    requests = [get_request(identifier) for identifier in identifiers]

    for r in tqdm(requests):
        r.result()

    for identifier in failed:
        tqdm.write(f"removing {identifier} from the dataset")
        for dataset in datasets:
            dataset.data = dataset.data.drop(
                dataset.data[
                    (dataset.data[PPI_NAME1] == identifier)
                    | (dataset.data[PPI_NAME2] == identifier)
                ].index
            )

    for dataset in datasets:
        print(f"{dataset.inpath} - Rows after filter sequence: {len(dataset.data)}")

        dataset.data.to_csv(dataset.outpath, index=False, sep=";")

    sequences_frame = pd.DataFrame.from_records(
        sequences, columns=["accession", "sequence"]
    )
    sequences_frame.to_csv(sequences_path, index=False, sep=";")


def ppi_trim(df):
    # filter on identifier columns
    df = df.filter(items=PPI_NAMES)

    print("Rows:", len(df))

    # remove rows with missing identifiers/wrong format identifiers
    df = df.loc[
        df[PPI_NAME1].str.startswith("uniprotkb:")
        & df[PPI_NAME2].str.startswith("uniprotkb:")
    ]
    df.drop_duplicates(PPI_NAMES, inplace=True)

    print("Rows after filter id and remove duplicates:", len(df))

    return df


@bacli.command
def ppi2(
    inpath1: str,
    outpath1: str = "PPI_positive.csv",
    inpath2: str = None,
    outpath2: str = "PPI_negative.csv",
    sequences_path: str = SEQUENCES_PATH,
):
    # Index        Protein_1_ID         protein_2_ID
    # 1 NP_663777.1  NP_001233.1
    # >NP_663777.1
    # MESSKKMDSPGALQTNPPLKLHTDRSAGTPVFVP...
    # >NP_001233.1
    # MARPHPWWLCVLGTLVGLSATPAPKSCPERHYWA...

    datasets = [(inpath1, outpath1)]
    if inpath2 is not None:
        datasets.append((inpath2, outpath2))

    sequences = dict()

    for inpath, outpath in datasets:
        interactions = list()
        to_remove = set()

        with open(inpath, "r") as file:
            header = file.readline()  # skip first line   # noqa: F841

            line = file.readline()
            while line:
                if line.startswith(">"):
                    id = line.lstrip(">").strip()  # noqa: A001
                    sequence = file.readline().strip()
                    # if sequence == "MLKSKTFLKKTRAGGVMKIVREHYLRDDIGCGAPGCAACGGAHEGPALEPQPQDPASSVCPQPHYLLPDTNVLLHQIDVLEDPAIRNVIVLQTVLQEVRNRSAPVYKRIRDVTNNQEKHFYTFTNEHHRETYVEQEQGENANDRNDRAIRVAAKWYNEHLKKMSADNQLQVIFITNDRRNKEKAIEEGIPAFTCEEYVKSLTANPELIDRLACLSEEGNEIESGKIIFSEHLPLSKLQQGIKSGTYLQGTFRASRENYLEATVWIHGDNEENKEIILQGLKHLNRAVHEDIVAVELLPKSQWVAPSSVVLHDEGQNEEDVEKEEETERMLKTAVSEKMLKPTGRVVGIIKRNWRPYCGMLSKSDIKESRRHLFTPADKRIPRIRIETRQASTLEGRRIIVAIDGWPRNSRYPNGHFVRNLGDVGEKETETEVLLLEHDVPHQPFSQAVLSFLPKMPWSITEKDMKNREDLRHLCICSVDPPGCTDIDDALHCRELENGNLEVGVHIADVSHFIRPGNALDQESARRGTTVYLCEKRIDMVPELLSSNLCSLKCDVDRLAFSCIWEMNHNAEILKTKFTKSVINSKASLTYAEAQLRIDSANMNDDITTSLRGLNKLAKILKKRRIEKGALTLSSPEVRFHMDSETHDPIDLQTKELRETNSMVEEFMLLANISVAKKIHEEFSEHALLRKHPAPPPSNYEILVKAARSRNLEIKTDTAKSLAESLDQAESPTFPYLNTLLRILATRCMMQAVYFCSGMDNDFHHYGLASPIYTHFTSPIRRYADVIVHRLLAVAIGADCTYPELTDKHKLADICKNLNFRHKMAQYAQRASVAFHTQLFFKSKGIVSEEAYILFVRKNAIVVLIPKYGLEGTVFFEEKDKPNPQLIYDDEIPSLKIEDTVFHVFDKVKVKIMLDSSNLQHQKIRMSLVEPQIPGISIPTDTSNMDLNGPKKKKMKLGK":
                    #     print(line)
                    #     print(id)
                    if id in sequences:
                        if sequence != sequences[id]:
                            print(
                                f"Sequences should match: \n{sequence}\n{sequences[id]}"
                            )
                        # assert sequence == sequences[id], f"Sequences should match: \n{sequence}\n{sequences[id]}"
                    else:
                        if len(sequence) > 6000:
                            if id not in to_remove:
                                print(f"skipping sequence with length {len(sequence)}")
                                to_remove.add(id)
                        else:
                            sequences[id] = sequence
                else:
                    index_id1, id2 = line.split("  ")
                    index, id1 = index_id1.split(" ")
                    interactions.append((id1.strip(), id2.strip()))

                line = file.readline()

        interactions = filter(
            lambda x: x[0] not in to_remove and x[1] not in to_remove, interactions
        )
        interactions_frame = pd.DataFrame.from_records(
            interactions, columns=[PPI_NAME1, PPI_NAME2]
        )
        interactions_frame.to_csv(outpath, index=False, sep=";")

    sequences_frame = pd.DataFrame.from_records(
        list(sequences.items()), columns=["identifier", "sequence"]
    )
    sequences_frame.to_csv(sequences_path, index=False, sep=";")
