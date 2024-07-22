"""
:author:    J.M. Algarín
:email:     josalggui@i3m.upv.es
:affiliation: MRILab, i3M, CSIC, Valencia, Spain

"""
import os
from datetime import datetime

from PyQt5.QtGui import QPixmap

from widgets.widget_toolbar_figures import FiguresToolBar


class FiguresController(FiguresToolBar):
    """
    Controller class for managing figures and screenshots.

    Args:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Inherits:
        FiguresToolBar: Base class for the figures toolbar.
    """
    def __init__(self, *args, **kwargs):
        """
        Initializes the FiguresController.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.
        """
        super(FiguresController, self).__init__(*args, **kwargs)

        self.action_full_screen.setCheckable(True)
        self.action_full_screen.triggered.connect(self.doFullScreen)
        self.action_screenshot.triggered.connect(self.doScreenshot)
        self.action_postprocessing.triggered.connect(self.openPostGui)

    def openPostGui(self):
        self.main.post_gui.showMaximized()
        self.main.post_gui.console.setup_console()

    def doFullScreen(self):
        """
        Toggles full-screen mode.

        Hides or shows specific GUI elements in full-screen mode.
        """
        if self.action_full_screen.isChecked():
            self.main.history_list.hide()
            self.main.sequence_list.hide()
            self.main.sequence_inputs.hide()
            self.main.console.hide()
            self.main.input_table.hide()
            self.main.custom_and_protocol.hide()
        else:
            self.main.history_list.show()
            self.main.sequence_list.show()
            self.main.sequence_inputs.show()
            self.main.console.show()
            self.main.input_table.show()
            self.main.custom_and_protocol.show()

    def doScreenshot(self):
        """
        Takes a screenshot of the main GUI and saves it.

        The screenshot is saved in the specified screenshot folder with a timestamp as the filename.
        """
        name = datetime.now()
        name_string = name.strftime("%Y.%m.%d.%H.%M.%S.%f")[:-3]
        file_name = name_string+".png"
        screenshot = QPixmap(self.main.size())
        self.main.render(screenshot)

        # Create screenshot folder
        screenshot_folder = self.main.session['directory'] + "/screenshots"
        if not os.path.exists(screenshot_folder):
            os.makedirs(screenshot_folder)

        # Save screenshot and print message
        screenshot.save(screenshot_folder+"/"+file_name)
        print("\nScreenshot saved in " + screenshot_folder+"/"+file_name)

class FiguresControllerPos(FiguresToolBar):
    """
    Controller class for managing figures and screenshots.

    Args:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Inherits:
        FiguresToolBar: Base class for the figures toolbar.
    """
    def __init__(self, *args, **kwargs):
        """
        Initializes the FiguresController.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.
        """
        super(FiguresControllerPos, self).__init__(*args, **kwargs)

        # Hide post button:
        self.action_postprocessing.setVisible(False)

        self.action_full_screen.setCheckable(True)
        self.action_full_screen.triggered.connect(self.doFullScreen)
        self.action_screenshot.triggered.connect(self.doScreenshot)

    def doFullScreen(self):
        """
        Toggles full-screen mode.

        Hides or shows specific GUI elements in full-screen mode.
        """
        if self.action_full_screen.isChecked():
            self.main.history_list.hide()
            self.main.methods_list.hide()
            self.main.tab_controller.hide()
            self.main.console.hide()
        else:
            self.main.history_list.show()
            self.main.methods_list.show()
            self.main.tab_controller.show()
            self.main.console.show()

    def doScreenshot(self):
        """
        Takes a screenshot of the main GUI and saves it.

        The screenshot is saved in the specified screenshot folder with a timestamp as the filename.
        """
        name = datetime.now()
        name_string = name.strftime("%Y.%m.%d.%H.%M.%S.%f")[:-3]
        file_name = name_string+".png"
        screenshot = QPixmap(self.main.size())
        self.main.render(screenshot)

        # Create screenshot folder
        screenshot_folder = self.main.session['directory'] + "/screenshots"
        if not os.path.exists(screenshot_folder):
            os.makedirs(screenshot_folder)

        # Save screenshot and print message
        screenshot.save(screenshot_folder+"/"+file_name)
        print("\nScreenshot saved in " + screenshot_folder+"/"+file_name)
