import json
from PySide6.QtWidgets import (QTableWidgetItem, QMessageBox, QApplication, QFileDialog, QComboBox, QLineEdit, QTableWidget, QRadioButton, QLabel, QCheckBox)
from Reload_Geometry import *

# ===============================================
# =========== SAVE/LOAD GUI PROGRESS ============
# ===============================================

def Save(self):
    # make user decide where to save file at
    self.save_file_name, _ = QFileDialog.getSaveFileName(None, "Save State", f'{self.directory}', "JSON File (*.json)")
    state = {"comboboxes": {},
             "line_edits": {},
             "tables": {},
             "radio_buttons": {},
             "labels": {},
             "checkboxes": {}}
             
    # change material definition of bulk forming solver to custom before saving
    self.INMaterialDefinition.setCurrentIndex(0)
             
    # Save widget states
    for attr_name, attr_value in self.__dict__.items():
        if isinstance(attr_value, QComboBox):
            state["comboboxes"][attr_name] = attr_value.currentText()
        elif isinstance(attr_value, QLineEdit):
            state["line_edits"][attr_name] = attr_value.text()
        elif isinstance(attr_value, QTableWidget):
            table_data = []
            for row in range(attr_value.rowCount()):
                row_data = []
                for col in range(attr_value.columnCount()):
                    item = attr_value.item(row, col)
                    row_data.append(item.text() if item else "")
                table_data.append(row_data)
            state["tables"][attr_name] = table_data
        elif isinstance(attr_value, QRadioButton):
            state["radio_buttons"][attr_name] = attr_value.isChecked()
        elif isinstance(attr_value, QLabel):
            state["labels"][attr_name] = attr_value.text()
        elif isinstance(attr_value, QCheckBox):
            state["checkboxes"][attr_name] = attr_value.isChecked()
            
    # Save the state to a file
    if self.save_file_name:
        with open(self.save_file_name, 'w') as file:
            json.dump(state, file)
            
def Load(self):
    # make user decide which file to load
    self.load_file_name, _ = QFileDialog.getOpenFileName(None, "Load State", f'{self.directory}', "JSON File (*.json)")
    
    if self.load_file_name:
        with open(self.load_file_name, 'r') as file:
            state = json.load(file)
        
        # Load widget states
        for attr_name, attr_value in self.__dict__.items():
            if isinstance(attr_value, QComboBox) and attr_name in state["comboboxes"]:
                index = attr_value.findText(state["comboboxes"][attr_name])
                if index >= 0:
                    attr_value.setCurrentIndex(index)
            elif isinstance(attr_value, QLineEdit) and attr_name in state["line_edits"]:
                attr_value.setText(state["line_edits"][attr_name])
            elif isinstance(attr_value, QTableWidget) and attr_name in state["tables"]:
                table_data = state["tables"][attr_name]
                attr_value.setRowCount(len(table_data))
                attr_value.setColumnCount(len(table_data[0]) if table_data else 0)
                for row, row_data in enumerate(table_data):
                    for col, cell_data in enumerate(row_data):
                        attr_value.setItem(row, col, QTableWidgetItem(cell_data))
            elif isinstance(attr_value, QRadioButton) and attr_name in state["radio_buttons"]:
                attr_value.setChecked(state["radio_buttons"][attr_name])
            elif isinstance(attr_value, QLabel) and attr_name in state["labels"]:
                attr_value.setText(state["labels"][attr_name])
            elif isinstance(attr_value, QCheckBox) and attr_name in state["checkboxes"]:
                attr_value.setChecked(state["checkboxes"][attr_name])
    
