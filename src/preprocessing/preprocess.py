""" Preprocess data files with pandas """
from itertools import chain

import pandas as pd
from requests_futures.sessions import (
    FuturesSession,
    ThreadPoolExecutor,
    # ProcessPoolExecutor,
)
from tqdm import tqdm

import src.bacli as bacli

# from collections import defaultdict
# from matplotlib import interactive


@bacli.command
def vdjdb(inpath: str, outpath: str):
    """ preprocess vdjdb file. Filter on TRB """

    data = pd.read_csv(inpath, sep="\t")

    # filter rows on TRB
    data = data.loc[data["gene"] == "TRB"]

    # filter columns
    columns = ["cdr3", "antigen.epitope"]
    data = data.filter(items=columns)

    data.to_csv(outpath, index=False, sep=";")


@bacli.command
def vdjdb_trb_human(inpath: str, outpath: str):
    """ preprocess vdjdb file. Filter on TRB """

    data = pd.read_csv(inpath, sep="\t")

    # filter rows on TRB
    data = data.loc[data["Gene"] == "TRB"]

    # filter rows on human
    data = data.loc[data["Species"] == "HomoSapiens"]

    # filter columns
    columns = ["CDR3", "Epitope"]
    data = data.filter(items=columns)

    # drop duplicates
    data = data.drop_duplicates(["CDR3", "Epitope"])

    data.to_csv(outpath, index=False, sep=";")


@bacli.command
def vdjdb_tra(inpath: str, outpath: str):
    """ preprocess vdjdb file. Filter on TRA """

    data = pd.read_csv(inpath, sep="\t")

    # filter rows on TRB
    data = data.loc[data["gene"] == "TRA"]

    # filter columns
    columns = ["cdr3", "antigen.epitope"]
    data = data.filter(items=columns)

    data.to_csv(outpath, index=False, sep=";")


@bacli.command
def vdjdb_hla(inpath: str, outpath: str):
    """ preprocess vdjdb file. Take a subset with mhc.a = HLA-A* """

    data = pd.read_csv(inpath, sep="\t")

    print(data.set_index(["gene", "mhc.a"]).count(level="mhc.a"))

    # filter rows on TRB
    data = data.loc[data["gene"] == "TRB"]
    data = data.loc[data["mhc.a"].str.startswith("HLA-A")]

    print(data.set_index(["gene", "mhc.a"]).count(level="mhc.a"))

    # filter columns
    columns = ["cdr3", "antigen.epitope"]
    data = data.filter(items=columns)

    data.to_csv(outpath, index=False, sep=";")


@bacli.command
def vdjdb_all(inpath: str, outpath: str):
    """ preprocess vdjdb file. No filtering """

    data = pd.read_csv(inpath, sep="\t")

    # filter columns
    columns = ["cdr3", "antigen.epitope"]
    data = data.filter(items=columns)

    data.to_csv(outpath, index=False, sep=";")


PPI_NAME1 = "Unique identifier for interactor A"
PPI_NAME2 = "Unique identifier for interactor B"
PPI_NAMES = [PPI_NAME1, PPI_NAME2]

SEQUENCES_PATH = "PPI_sequences.csv"


@bacli.command
def ppi(
    inpath1: str,
    outpath1: str,
    inpath2: str = None,
    outpath2: str = None,
    sequencesPath: str = SEQUENCES_PATH,
):
    def getAccession(identifier):
        if not identifier.startswith("uniprotkb:"):
            return None
        accession = identifier.split(":")[1]
        # DO NOT DROP PART AFTER dash (-)
        # accession = accession.split('-')[0]
        return accession

    def getUrl(identifier):
        accession = getAccession(identifier)
        url = f"https://www.ebi.ac.uk/proteins/api/proteins/{accession}"
        # print(url)
        return url

    class PpiDataset(object):
        def __init__(self, inpath, outpath):
            self.inpath = inpath
            self.outpath = outpath
            self.data = ppi_trim(pd.read_csv(self.inpath, sep="\t"))

        def getIdentifiers(self):
            return set(chain(self.data[PPI_NAME1], self.data[PPI_NAME2]))

    datasets = [PpiDataset(inpath1, outpath1)]
    if inpath2 is not None:
        datasets.append(PpiDataset(inpath2, outpath2))

    identifiers = set(
        identifier for dataset in datasets for identifier in dataset.getIdentifiers()
    )
    print("Unique identifiers: ", len(identifiers))

    sequences = list()
    failed = list()

    def getResponseHandler(identifier):
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
            jsonData = response.json()
            sequence = jsonData["sequence"]["sequence"]

            sequences.append((identifier, sequence))

        return handler

    session = FuturesSession(executor=ThreadPoolExecutor(max_workers=20))

    def getRequest(identifier):
        url = getUrl(identifier)
        return session.get(url, hooks={"response": getResponseHandler(identifier)})

    requests = [getRequest(identifier) for identifier in identifiers]

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

    sequencesFrame = pd.DataFrame.from_records(
        sequences, columns=["accession", "sequence"]
    )
    sequencesFrame.to_csv(sequencesPath, index=False, sep=";")


def ppi_trim(dataFrame):
    # filter on identifier columns
    dataFrame = dataFrame.filter(items=PPI_NAMES)

    print("Rows:", len(dataFrame))

    # remove rows with missing identifiers/wrong format identifiers
    dataFrame = dataFrame.loc[
        dataFrame[PPI_NAME1].str.startswith("uniprotkb:")
        & dataFrame[PPI_NAME2].str.startswith("uniprotkb:")
    ]
    dataFrame.drop_duplicates(PPI_NAMES, inplace=True)

    print("Rows after filter id and remove duplicates:", len(dataFrame))

    return dataFrame


@bacli.command
def ppi2(
    inpath1: str,
    outpath1: str = "PPI_positive.csv",
    inpath2: str = None,
    outpath2: str = "PPI_negative.csv",
    sequencesPath: str = SEQUENCES_PATH,
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
        toRemove = set()

        with open(inpath, "r") as file:
            header = file.readline()

            line = file.readline()
            while line:
                if line.startswith(">"):
                    id = line.lstrip(">").strip()
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
                            if id not in toRemove:
                                print(f"skipping sequence with length {len(sequence)}")
                                toRemove.add(id)
                        else:
                            sequences[id] = sequence
                else:
                    index_id1, id2 = line.split("  ")
                    index, id1 = index_id1.split(" ")
                    interactions.append((id1.strip(), id2.strip()))

                line = file.readline()

        interactions = filter(
            lambda x: x[0] not in toRemove and x[1] not in toRemove, interactions
        )
        interactionsFrame = DataFrame.from_records(
            interactions, columns=[PPI_NAME1, PPI_NAME2]
        )
        interactionsFrame.to_csv(outpath, index=False, sep=";")

    sequencesFrame = DataFrame.from_records(
        list(sequences.items()), columns=["identifier", "sequence"]
    )
    sequencesFrame.to_csv(sequencesPath, index=False, sep=";")
