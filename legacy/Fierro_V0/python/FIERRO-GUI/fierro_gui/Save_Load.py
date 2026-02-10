import json
from PySide6.QtWidgets import (QTableWidgetItem, QMessageBox, QApplication, QFileDialog, QComboBox, QLineEdit, QTableWidget, QRadioButton, QLabel, QCheckBox, QTreeWidget, QTreeWidgetItem, QStackedWidget, QWidget, QPushButton)
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
             "checkboxes": {},
             "tree_widgets": {},
             "stacked_widgets": {}}
             
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
        elif isinstance(attr_value, QTreeWidget):
            tree_data = []
            for i in range(attr_value.topLevelItemCount()):
                item = attr_value.topLevelItem(i)
                tree_data.append(get_tree_widget_items(item))
            state["tree_widgets"][attr_name] = tree_data
        elif isinstance(attr_value, QStackedWidget):
            stacked_widget_data = []
            for i in range(attr_value.count()):
                page = attr_value.widget(i)
                page_data = {
                    "type": type(page).__name__,
                    "widgets": {}
                }
                for child in page.findChildren(QWidget):
                    if isinstance(child, QComboBox):
                        page_data["widgets"][child.objectName()] = {"type": "QComboBox", "value": child.currentText()}
                    elif isinstance(child, QLineEdit):
                        page_data["widgets"][child.objectName()] = {"type": "QLineEdit", "value": child.text()}
                    elif isinstance(child, QLabel):  # Save QLabel
                        page_data["widgets"][child.objectName()] = {
                            "type": "QLabel",
                            "text": child.text()
                        }
                    elif isinstance(child, QTableWidget):  # Save QTableWidget
                        table_data = []
                        for row in range(child.rowCount()):
                            row_data = []
                            for col in range(child.columnCount()):
                                item = child.item(row, col)
                                row_data.append(item.text() if item else "")
                            table_data.append(row_data)
                        page_data["widgets"][child.objectName()] = {
                            "type": "QTableWidget",
                            "data": table_data
                        }
                    elif isinstance(child, QPushButton):
                        page_data["widgets"][child.objectName()] = {
                            "type": "QPushButton",
                            "enabled": child.isEnabled()
                        }
                    # ... add other widget types as needed
                stacked_widget_data.append(page_data)
            state["stacked_widgets"][attr_name] = stacked_widget_data

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
            elif isinstance(attr_value, QTreeWidget) and attr_name in state["tree_widgets"]:
                attr_value.clear()
                load_tree_widget_items(attr_value, state["tree_widgets"][attr_name])
            elif isinstance(attr_value, QStackedWidget) and attr_name in state["stacked_widgets"]:
                stacked_widget_data = state["stacked_widgets"][attr_name]
                for page_data in stacked_widget_data:
                    page = globals()[page_data["type"]]()  # Create new page of the saved type
                    for widget_name, widget_data in page_data["widgets"].items():
                        widget = getattr(page, widget_name, None)
                        if widget:
                            if widget_data["type"] == "QComboBox":
                                index = widget.findText(widget_data["value"])
                                if index >= 0:
                                    widget.setCurrentIndex(index)
                            elif widget_data["type"] == "QLineEdit":
                                widget.setText(widget_data["value"])
                            elif widget_data["type"] == "QLabel":  # Load QLabel
                                widget.setText(widget_data["text"])
                            elif widget_data["type"] == "QTableWidget":  # Load QTableWidget
                                table = widget  # Assuming 'widget' is a QTableWidget
                                table.setRowCount(len(widget_data["data"]))
                                table.setColumnCount(len(widget_data["data"][0]) if widget_data["data"] else 0)
                                
                                for row_index, row_data in enumerate(widget_data["data"]):
                                    for col_index, cell_data in enumerate(row_data):
                                        table.setItem(row_index, col_index, QTableWidgetItem(cell_data))
                            elif widget_data["type"] == "QPushButton":
                                widget.setEnabled(widget_data["enabled"])
                    # ... handle other widget types
                    attr_value.addWidget(page)

def get_tree_widget_items(item):
    data = {
        "text": [item.text(i) for i in range(item.columnCount())],
        "children": []
    }
    for i in range(item.childCount()):
        data["children"].append(get_tree_widget_items(item.child(i)))
    return data
        
def load_tree_widget_items(tree_widget, items_data):
    for item_data in items_data:
        item = QTreeWidgetItem(tree_widget)
        set_tree_widget_item_data(item, item_data)

def set_tree_widget_item_data(item, data):
    for i, text in enumerate(data["text"]):
        item.setText(i, text)
    for child_data in data["children"]:
        child_item = QTreeWidgetItem(item)
        set_tree_widget_item_data(child_item, child_data)
    
