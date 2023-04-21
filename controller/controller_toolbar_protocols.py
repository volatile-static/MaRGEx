"""
@author:    José Miguel Algarín
@email:     josalggui@i3m.upv.es
@affiliation:MRILab, i3M, CSIC, Valencia, Spain
"""
import csv
import os
import platform
import shutil

from PyQt5.QtWidgets import QFileDialog

from seq.sequences import defaultsequences
from widgets.widget_toolbar_protocols import ProtocolsToolBar


class ProtocolsController(ProtocolsToolBar):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.action_new_protocol.triggered.connect(self.newProtocol)
        self.action_del_protocol.triggered.connect(self.delProtocol)
        self.action_new_sequence.triggered.connect(self.newSequence)
        self.action_del_sequence.triggered.connect(self.delSequence)

    def delProtocol(self):
        # Open a file dialog to get the filename to save to
        directory = 'protocols'
        folder_name = QFileDialog.getExistingDirectory(self.main, "Remove protocol", directory)

        if folder_name:
            shutil.rmtree(folder_name)
            print("\nProtocol removed")
            self.main.protocol_list.updateProtocolList()

    def delSequence(self):
        # Get the current protocol
        protocol = self.main.protocol_list.getCurrentProtocol()

        # Open a file dialog to get the filename to save to
        directory = 'protocols/%s' % protocol
        file_name, _ = QFileDialog.getOpenFileName(self.main, 'Remove sequence from protocol', directory, '(*.csv)')

        # Delete protocol
        if file_name:
            os.remove(file_name)
            print("Protocol removed")
            self.main.protocol_inputs.updateProtocolInputs()

    def newProtocol(self):
        # Open a file dialog to get the filename to save to
        file_name, _ = QFileDialog.getSaveFileName(self.main, 'New Protocol', 'protocols', '')

        if file_name:
            # Delete extension
            file_name = file_name.split('.')[0]

            # Check if the folder is the good one
            directory = os.path.dirname(file_name).split('/')[-1]
            if directory != 'protocols':
                print("Error. New protocols should be in 'protocols' folder.")
                return

            if not os.path.exists(file_name):
                os.makedirs(file_name)
                print("New protocol created successfully")
                self.main.protocol_list.updateProtocolList()
            else:
                print("Protocol already exist")

    def newSequence(self):
        # Get the current protocol
        protocol = self.main.protocol_list.getCurrentProtocol()

        # Get the current sequence
        seq_name = self.main.sequence_list.getCurrentSequence()
        sequence = defaultsequences[seq_name]

        # Open a file dialog to get the filename to save to
        directory = 'protocols/%s' % protocol
        file_name, _ = QFileDialog.getSaveFileName(self.main, 'Add sequence to protocol', directory, '(*.csv)')
        if file_name:
            if platform.system()=='Linux':
                file_name = "%s_%s.csv" % (seq_name, file_name.split('/')[-1])
            else:
                file_name = "%s_%s" % (seq_name, file_name.split('/')[-1])

            # Save csv with input parameters
            with open('%s/%s' % (directory, file_name), 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=sequence.mapKeys)
                writer.writeheader()
                map_vals = {}
                for key in sequence.mapKeys:  # take only the inputs from mapVals
                    map_vals[key] = sequence.mapVals[key]
                writer.writerows([sequence.mapNmspc, map_vals])

            self.main.protocol_inputs.updateProtocolInputs()

            print("\n%s sequence added to the %s protocol" % (file_name.split('.')[0], protocol))
