import csv

import pandas as pd

CAT = ['1', '2', '3', '4']


def count_cat(production_matrix: str):
    """Count occurrences of compounds for each type of production

    Args:
        production_matrix (str): tsv file indicating for each holobiont (columns) and each compound (rows) which way
        the compound is produced.

    Returns:
        dict: dictionary
    """
    df = pd.DataFrame(columns=['Compound'] + CAT)
    with open(production_matrix, 'r') as fi:
        lines = csv.reader(fi, delimiter='\t')
        lines.__next__()
        for l in lines:
            row = [l[0]] + [l.count(c) for c in CAT]
            df2 = pd.DataFrame([row], columns=['Compound'] + CAT)
            df = pd.concat([df, df2])
    return df[0:20]

    # df = pd.DataFrame(columns=['Compound', 'Cat', 'Count'])
    # with open(production_matrix, 'r') as fi:
    #     lines = csv.reader(fi, delimiter='\t')
    #     lines.__next__()
    #     for l in lines:
    #         for c in CAT:
    #             row = [l[0], c, l.count(c)]
    #             df2 = pd.DataFrame([row], columns=['Compound', 'Cat', 'Count'])
    #             df = pd.concat([df, df2])
    # return df

    # cat_count = dict()
    # with open(production_matrix, 'r') as fi:
    #     lines = csv.reader(fi, delimiter='\t')
    #     lines.__next__()
    #     for l in lines:
    #         cat_count[l[0]] = dict()
    #         for cat in CAT:
    #             cat_count[l[0]][cat] = l.count(cat)
    # print(cat_count)
    # return cat_count


def bar_occurrence_cat(cat_count):
    import plotly.express as px
    df = px.data.tips()
    print(df)
    # fig = px.bar(df, x="sex", y="total_bill", color='time')
    # fig.show()

    print(cat_count)
    fig = px.bar(cat_count, x='2', y='Compound', orientation='h')
    fig.show()


mat_file = '../../example/outputs/heatmap/solo/run_matrix.tsv'
d = count_cat(mat_file)
bar_occurrence_cat(d)