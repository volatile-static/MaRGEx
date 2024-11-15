"""
:author:    J.M. Algarín
:email:     josalggui@i3m.upv.es
:affiliation: MRILab, i3M, CSIC, Valencia, Spain

"""
import subprocess
import sys
import threading

from PyQt5.QtCore import QEvent

from seq.sequences import defaultsequences
from ui.window_main import MainWindow


class MainController(MainWindow):
    def __init__(self, *args, **kwargs):
        super(MainController, self).__init__(*args, **kwargs)

        self.set_session(self.session)

        self.initializeThread()

    def set_demo(self, demo):
        self.demo = demo

    def set_session(self, session):
        # Get repo version
        try:
            tag = subprocess.check_output(['git', 'describe', '--tags'], stderr=subprocess.STDOUT).strip().decode(
                'utf-8')
        except subprocess.CalledProcessError as e:
            print(f"Error getting Git tag: {e.output.decode('utf-8')}")
            tag = ""

        # Set window title
        self.session = session
        self.setWindowTitle("MaRGE " + tag + ": " + session['directory'])
        # Add the session to all sequences
        for sequence in defaultsequences.values():
            sequence.session = session

    def initializeThread(self):
        # Start the sniffer
        thread = threading.Thread(target=self.history_list.waitingForRun)
        thread.start()
        print("Sniffer initialized.\n")

    def set_console(self):
        self.layout_left.addWidget(self.console)

    def changeEvent(self, event):
        if event.type() == QEvent.ActivationChange:  # Event type 99
            if self.isActiveWindow():
                self.set_console()
        super().changeEvent(event)

    def closeEvent(self, event):
        """
        Shuts down the application on close.

        This method is called when the application is being closed. It sets the `app_open` flag to False, restores
        `sys.stdout` to its default value, and performs additional cleanup tasks if the `demo` flag is not set.
        It also prints a closing message to the console.

        Args:
            event (QCloseEvent): The close event triggered by the user.

        Returns:
            None
        """
        # Return stdout to defaults.
        sys.stdout = sys.__stdout__
            
        print('\nMain GUI closed successfully!')

        self.parent.show()

        super().closeEvent(event)
