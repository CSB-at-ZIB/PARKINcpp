'''
Some helper functions to handle a data file.
'''
import csv
import logging
from stabledict import StableDict
from parkin import MeasurementList, MeasurementPoint

__author__ = 'Moritz Wade'

def _register_csv_dialect():
    csv.register_dialect("ExperimentalData", delimiter=';', quotechar='"', skipinitialspace = True)

def save_to_file(filename, timepoints, species_ids, species_data):
    _register_csv_dialect()

    num_species = len(species_ids)
    num_entries = len(timepoints)

    try:
        writer = csv.writer(open(filename, "w"), dialect="ExperimentalData")

        writer.writerow(["Timepoint"].extend(species_ids))

        for row in xrange(num_entries):
            rowData = []
            rowData.append(timepoints[row])   # adds the timepoint up front
            for col in xrange(num_species):
                index = row * 2 + col
                species_value = species_data[index]
                rowData.append(species_value)
            writer.writerow(rowData)

    except Exception, e:
        logging.error("Error while trying to write CSV file: %s\nError: %s" % (filename,e))
    pass

def load_file(filename):
    _register_csv_dialect()
    reader = csv.reader(open(filename), dialect="SimulationData")

    timepoints = []
    measurement_list = []


    for (i, row) in enumerate(reader):  # we don't use the i... yet
        values = StableDict()
        for (j,value) in enumerate(row):
            if j == 0:  # jump over the 0th row
                time = value # time value is always in the 0th column
                continue
            species_id = reader.fieldnames[j]
            values[species_id] = value
        timepoints.append(time)
        measurement_list.append(values)

    measurement_vector = MeasurementList(len(measurement_list))
    for i, (time, values) in enumerate(measurement_list):
        meas_point = MeasurementPoint()
        for id, value in values.items():
            meas_point[id] = (value, 0.0)   # dummy weight of 0.0 for now
        measurement_vector[i] = meas_point


    return (timepoints, measurement_vector)
