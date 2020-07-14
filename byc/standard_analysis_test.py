import pandas
from byc.standard_analysis import *

def test_Cell_Data():
    pass


def test_set_processed_traces():
    window_width = 3
    collection_interval = 10
    processed_traces, trace_filenames = set_processed_traces(window_width,
    														 collection_interval,
    														 manual_select=False)

    trace_types = [type(trace) for trace in processed_traces]
    type_matches = [datatype == pandas.core.frame.DataFrame for datatype in trace_types]

    assert False not in type_matches, "standard_analysis.set_processed_traces is returning non DataFrame trace types"


