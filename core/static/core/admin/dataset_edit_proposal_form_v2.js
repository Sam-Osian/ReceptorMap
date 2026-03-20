(function () {
    var NEW_CHOICE_VALUE = "__new__";

    function bySelector(selector) {
        return document.querySelector(selector);
    }

    function field(name) {
        return document.getElementById("id_" + name);
    }

    function rowForField(name) {
        return bySelector(".form-row.field-" + name);
    }

    function asText(value) {
        if (value === null || value === undefined) {
            return "";
        }
        return String(value);
    }

    function setValue(name, value) {
        var el = field(name);
        if (!el) {
            return;
        }
        if (el.type === "checkbox") {
            el.checked = Boolean(value);
            return;
        }
        el.value = asText(value);
    }

    function setEnabled(name, enabled) {
        var el = field(name);
        if (el) {
            el.disabled = !enabled;
        }
    }

    function setVisible(el, visible) {
        if (el) {
            el.style.display = visible ? "" : "none";
        }
    }

    function setRowVisible(name, visible) {
        setVisible(rowForField(name), visible);
    }

    function choiceHasOption(selectEl, value) {
        if (!selectEl) {
            return false;
        }
        for (var i = 0; i < selectEl.options.length; i += 1) {
            if (selectEl.options[i].value === value) {
                return true;
            }
        }
        return false;
    }

    function setSelectOrNew(selectName, newName, value) {
        var selectEl = field(selectName);
        var newEl = field(newName);
        var text = asText(value);
        if (!selectEl || !newEl) {
            return;
        }

        if (text && choiceHasOption(selectEl, text)) {
            selectEl.value = text;
            newEl.value = "";
        } else if (text) {
            selectEl.value = NEW_CHOICE_VALUE;
            newEl.value = text;
        } else {
            selectEl.value = "";
            newEl.value = "";
        }
    }

    function syncNewFieldVisibility(selectName, newName) {
        var selectEl = field(selectName);
        var newEl = field(newName);
        if (!selectEl || !newEl) {
            return;
        }
        var show = selectEl.value === NEW_CHOICE_VALUE;
        setRowVisible(newName, show);
        setEnabled(newName, show);
        if (!show) {
            newEl.value = "";
        }
    }

    function snapshotText(row) {
        if (!row) {
            return "";
        }
        var columns = [
            "drug",
            "drug_class",
            "target",
            "activity",
            "activity_recoded",
            "ki_lower_nm",
            "ki_upper_nm",
            "modelled_ki_nm",
            "ki_is_range",
            "log10_ki_nm",
            "pKi",
        ];
        return columns
            .map(function (column) {
                return column + ": " + asText(row[column]);
            })
            .join("\n");
    }

    document.addEventListener("DOMContentLoaded", function () {
        var rowSelect = field("existing_row");
        var actionSelect = field("action");
        var snapshot = field("current_row_snapshot");
        if (!rowSelect) {
            return;
        }

        var rowSelectionSection = bySelector(".workflow-row-selection");
        var rowValuesSection = bySelector(".workflow-row-values");
        var evidenceSection = bySelector(".workflow-evidence");
        var researchSection = bySelector(".workflow-research-context");
        var deleteSection = bySelector(".workflow-delete-details");

        var rowMap = {};
        try {
            rowMap = JSON.parse(rowSelect.dataset.rowMap || "{}");
        } catch (err) {
            rowMap = {};
        }
        var saveButtons = Array.from(
            document.querySelectorAll(
                "input[name='_save'], input[name='_addanother'], input[name='_continue']"
            )
        );
        var duplicateWarning = document.createElement("div");
        duplicateWarning.className = "help";
        duplicateWarning.style.color = "#b91c1c";
        duplicateWarning.style.marginTop = "0.35rem";
        duplicateWarning.style.display = "none";
        duplicateWarning.textContent =
            "This drug/target combination already exists. Use the Update workflow instead.";
        if (rowValuesSection) {
            rowValuesSection.insertBefore(duplicateWarning, rowValuesSection.firstChild);
        } else {
            var anchorRow = rowForField("target_new") || rowForField("target_select");
            if (anchorRow && anchorRow.parentNode) {
                anchorRow.parentNode.insertBefore(duplicateWarning, anchorRow.nextSibling);
            }
        }
        var updateSummary = document.createElement("div");
        updateSummary.className = "help";
        updateSummary.style.marginTop = "0.65rem";
        updateSummary.style.padding = "0.55rem 0.65rem";
        updateSummary.style.border = "1px solid #d7deef";
        updateSummary.style.borderRadius = "0.45rem";
        updateSummary.style.background = "#f8faff";
        updateSummary.style.display = "none";
        if (rowValuesSection) {
            rowValuesSection.appendChild(updateSummary);
        }

        function currentAddValue(selectName, newName) {
            var selectEl = field(selectName);
            var newEl = field(newName);
            if (!selectEl) {
                return "";
            }
            if (selectEl.value === NEW_CHOICE_VALUE) {
                return newEl ? (newEl.value || "").trim() : "";
            }
            return (selectEl.value || "").trim();
        }

        function addComboExists() {
            var drugValue = currentAddValue("drug_select", "drug_new");
            var targetValue = currentAddValue("target_select", "target_new");
            if (!drugValue || !targetValue) {
                return false;
            }
            return Object.keys(rowMap).some(function (key) {
                var row = rowMap[key];
                return row && row.drug === drugValue && row.target === targetValue;
            });
        }

        function setSaveEnabled(enabled) {
            saveButtons.forEach(function (btn) {
                btn.disabled = !enabled;
            });
        }

        function refreshDuplicateState() {
            var isAdd = currentMode() === "add";
            var duplicate = isAdd && addComboExists();
            duplicateWarning.style.display = duplicate ? "" : "none";
            setSaveEnabled(!duplicate);
        }

        function getSelectOrNewValue(selectName, newName) {
            var selectEl = field(selectName);
            var newEl = field(newName);
            if (!selectEl) {
                return "";
            }
            if (selectEl.value === NEW_CHOICE_VALUE) {
                return newEl ? (newEl.value || "").trim() : "";
            }
            return (selectEl.value || "").trim();
        }

        function valuesDiffer(oldValue, newValue, numeric) {
            if (numeric) {
                var oldNum = parseFloat(asText(oldValue));
                var newNum = parseFloat(asText(newValue));
                if (Number.isNaN(oldNum) && Number.isNaN(newNum)) {
                    return false;
                }
                if (!Number.isNaN(oldNum) && !Number.isNaN(newNum)) {
                    return Math.abs(oldNum - newNum) > 1e-9;
                }
                return true;
            }
            return asText(oldValue).trim() !== asText(newValue).trim();
        }

        function escapeHtml(text) {
            return asText(text)
                .replace(/&/g, "&amp;")
                .replace(/</g, "&lt;")
                .replace(/>/g, "&gt;")
                .replace(/\"/g, "&quot;")
                .replace(/'/g, "&#39;");
        }

        function refreshUpdateSummary() {
            var isUpdate = currentMode() === "update";
            var key = rowSelect.value || "";
            var row = rowMap[key];
            if (!isUpdate || !row) {
                updateSummary.style.display = "none";
                updateSummary.innerHTML = "";
                return;
            }

            var changes = [];
            var newDrugClass = getSelectOrNewValue("drug_class_select", "drug_class_new");
            var newActivity = getSelectOrNewValue("activity_select", "activity_new");
            var newKiLower = (field("ki_lower_nm") ? field("ki_lower_nm").value : "").trim();
            var newKiUpper = (field("ki_upper_nm") ? field("ki_upper_nm").value : "").trim();

            if (valuesDiffer(row.drug_class, newDrugClass, false)) {
                changes.push("Drug class: " + asText(row.drug_class) + " -> " + newDrugClass);
            }
            if (valuesDiffer(row.activity, newActivity, false)) {
                changes.push("Activity: " + asText(row.activity) + " -> " + newActivity);
            }
            if (valuesDiffer(row.ki_lower_nm, newKiLower, true)) {
                changes.push("Ki lower bound (nM): " + asText(row.ki_lower_nm) + " -> " + newKiLower);
            }
            if (valuesDiffer(row.ki_upper_nm, newKiUpper, true)) {
                changes.push("Ki upper bound (nM): " + asText(row.ki_upper_nm) + " -> " + newKiUpper);
            }

            updateSummary.style.display = "";
            if (!changes.length) {
                updateSummary.innerHTML = "<strong>Changed fields:</strong> None yet.";
                return;
            }
            updateSummary.innerHTML =
                "<strong>Changed fields:</strong><ul style='margin:0.35rem 0 0 1.1rem;'>" +
                changes
                    .map(function (item) {
                        return "<li>" + escapeHtml(item) + "</li>";
                    })
                    .join("") +
                "</ul>";
        }

        function currentMode() {
            return actionSelect ? actionSelect.value : "";
        }

        function applyRowSelection(forceFillValues) {
            var key = rowSelect.value || "";
            var row = rowMap[key];
            if (!row) {
                if (snapshot) {
                    snapshot.value = "";
                }
                return;
            }

            if (snapshot) {
                snapshot.value = snapshotText(row);
            }

            setSelectOrNew("drug_select", "drug_new", row.drug);
            setSelectOrNew("target_select", "target_new", row.target);
            if (forceFillValues) {
                setSelectOrNew("drug_class_select", "drug_class_new", row.drug_class);
                setSelectOrNew("activity_select", "activity_new", row.activity);
                setValue("ki_lower_nm", row.ki_lower_nm);
                setValue("ki_upper_nm", row.ki_upper_nm);
            }
        }

        function setModeState() {
            var action = currentMode();
            var isAdd = action === "add";
            var isUpdate = action === "update";
            var isDelete = action === "delete";

            setVisible(rowSelectionSection, isUpdate || isDelete);
            setVisible(rowValuesSection, isAdd || isUpdate);
            setVisible(evidenceSection, isAdd || isUpdate || isDelete);
            setVisible(researchSection, isAdd || isUpdate);
            setVisible(deleteSection, isDelete);
            setRowVisible("current_row_snapshot", isDelete);

            setEnabled("existing_row", isUpdate || isDelete);
            setEnabled("current_row_snapshot", isUpdate || isDelete);

            // Drug/target are only user-selectable during Add.
            setRowVisible("drug_select", isAdd);
            setRowVisible("drug_new", isAdd);
            setRowVisible("target_select", isAdd);
            setRowVisible("target_new", isAdd);
            setEnabled("drug_select", isAdd);
            setEnabled("target_select", isAdd);

            // Shared Add/Update editable fields.
            setEnabled("drug_class_select", isAdd || isUpdate);
            setEnabled("activity_select", isAdd || isUpdate);
            setEnabled("ki_lower_nm", isAdd || isUpdate);
            setEnabled("ki_upper_nm", isAdd || isUpdate);

            if (isAdd) {
                rowSelect.value = "";
                if (snapshot) {
                    snapshot.value = "";
                }
            }

            if ((isUpdate || isDelete) && rowSelect.value) {
                applyRowSelection(isUpdate);
            }

            syncNewFieldVisibility("drug_select", "drug_new");
            syncNewFieldVisibility("target_select", "target_new");
            syncNewFieldVisibility("drug_class_select", "drug_class_new");
            syncNewFieldVisibility("activity_select", "activity_new");
            refreshDuplicateState();
            refreshUpdateSummary();
        }

        rowSelect.addEventListener("change", function () {
            applyRowSelection(true);
            setModeState();
        });

        ["drug_select", "target_select", "drug_class_select", "activity_select"].forEach(function (name) {
            var selectEl = field(name);
            if (!selectEl) {
                return;
            }
            selectEl.addEventListener("change", function () {
                syncNewFieldVisibility(name, name.replace("_select", "_new"));
                refreshDuplicateState();
                refreshUpdateSummary();
            });
        });
        ["drug_new", "target_new", "drug_class_new", "activity_new", "ki_lower_nm", "ki_upper_nm"].forEach(function (name) {
            var inputEl = field(name);
            if (!inputEl) {
                return;
            }
            inputEl.addEventListener("input", function () {
                refreshDuplicateState();
                refreshUpdateSummary();
            });
        });

        if (actionSelect) {
            actionSelect.addEventListener("change", setModeState);
        }

        setModeState();
    });
})();
