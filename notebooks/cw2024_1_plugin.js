/**
 * @file cw2024_1_plugin.js
 * @author Christophe Avenel
 */

/**
 * @namespace cw2024_1_plugin
 * @classdesc The root namespace for cw2024_1_plugin.
 */
var cw2024_1_plugin;
cw2024_1_plugin = {
  name: "Spatialdata export Plugin",
  parameters: {
    _labels: {
      label: "Labels in JSON format",
      type: "text",
      default: "[]"
    },
    _addLayers: {
      label: "Add layer names",
      type: "checkbox",
      default: true
    },
    _layerFontSize: {
      label: "Layer name font size",
      type: "number",
      default: 12
    },
    _fontColor: {
      label: "Font color",
      type: "text",
      default: "#000000"
    },
  },
};

/**
 * @summary */
cw2024_1_plugin.init = async function (container) {
  await overlayUtils.waitLayersReady();
  cw2024_1_plugin.refresh();
  tmapp["ISS_viewer"].addHandler("animation", function animationFinishHandler(event){
    let classOverlays = document.getElementsByClassName("cw2024_1_plugin_overlay");
    for (let i = classOverlays.length - 1; i >= 0; i--) {
      // get zoomable attribute
      let zoomable = classOverlays[i].getAttribute("data-zoomable");
      if (zoomable == "true") {
        // get inner div of class cw2024_1_plugin_div
        let cw2024_1_plugin_div = classOverlays[i].getElementsByClassName("cw2024_1_plugin_div")[0];
        cw2024_1_plugin_div.style.zoom = tmapp.ISS_viewer.viewport.getZoom(false);
      }
    }
  });
};

cw2024_1_plugin.inputTrigger = function (parameterName) {
  cw2024_1_plugin.refresh();
}

cw2024_1_plugin.addOverlay = async function (x,y,text, size, rotation, zoomable, placement) {
    const OSDviewer = tmapp.ISS_viewer;
    let boundBoxOverlay = document.createElement("div");
    let rand_id = Math.random().toString(36).substring(7);
    boundBoxOverlay.className = "cw2024_1_plugin_overlay";
    boundBoxOverlay.style.fontSize = size + "px";
    boundBoxOverlay.style.zIndex = 1000;
    boundBoxOverlay.style.color = cw2024_1_plugin.get("_fontColor");
    boundBoxOverlay.style.width = "0px";
    boundBoxOverlay.style.whiteSpace = "nowrap";
    // add data-zoomable attribute
    boundBoxOverlay.setAttribute("data-zoomable", zoomable);
    boundBoxOverlay.innerHTML = `<div class='cw2024_1_plugin_div' style='transform:rotate(${rotation}deg)' id='cw2024_1_plugin_div_${rand_id}'>${text}</div>`;
    // Add random ID
    boundBoxOverlay.id = "cw2024_1_plugin_overlay_" + Math.random().toString(36).substring(7);
    let position = new OpenSeadragon.Point(x, y);
    OSDviewer.addOverlay(boundBoxOverlay, position, placement);
    if (zoomable) {
        // get inner div of class cw2024_1_plugin_div
        let cw2024_1_plugin_div = boundBoxOverlay.getElementsByClassName("cw2024_1_plugin_div")[0];
        cw2024_1_plugin_div.style.zoom = tmapp.ISS_viewer.viewport.getZoom(false);
    }
}

cw2024_1_plugin.refresh = function () {
    const OSDviewer = tmapp.ISS_viewer;
    // Remove all overlays of class cw2024_1_plugin_overlay that are not used
    let classOverlays = document.getElementsByClassName("cw2024_1_plugin_overlay");
    for (let i = classOverlays.length - 1; i >= 0; i--) {
        console.log(i);
        OSDviewer.removeOverlay(classOverlays[i])
    }
    if (cw2024_1_plugin.get("_addLayers")) {
        var count = OSDviewer.world.getItemCount();
        for (var i = 0; i < count; i++) {
            let item = OSDviewer.world.getItemAt(i);
            let layer = tmapp.layers[i];
            if (item.source) {
                let position = item.imageToViewportCoordinates(0,0)
                cw2024_1_plugin.addOverlay(position.x, position.y, layer.name, cw2024_1_plugin.get("_layerFontSize"), 0, false, OpenSeadragon.Placement.BOTTOM_LEFT);
            }
        }
    }
    JSON.parse(cw2024_1_plugin.get("_labels")).forEach(function (label) {
        cw2024_1_plugin.addOverlay(label.x, label.y, label.text, label.size, label.rotation, label.zoomable, OpenSeadragon.Placement[label.placement]);
    });
    tmapp.ISS_viewer.viewport.zoomBy(1)


};

/*

Example of json format for labels:
[
    {
        "x": 0.5,
        "y": 0.5,
        "text": "Center",
        "size": 20,
        "rotation": 90,
        "zoomable": true,
        "placement": "CENTER"
    },
    {
        "x": 0.1,
        "y": 0.1,
        "text": "Top left",
        "size": 20,
        "rotation": 0,
        "zoomable": true,
        "placement": "TOP_LEFT"
    },
    {
        "x": 0.9,
        "y": 0.1,
        "text": "Top right",
        "size": 20,
        "rotation": 0,
        "zoomable": true,
        "placement": "TOP_RIGHT"
    },
    {
        "x": 0.1,
        "y": 0.9,
        "text": "Bottom left",
        "size": 20,
        "rotation": 0,
        "zoomable": true,
        "placement": "BOTTOM_LEFT"
    },
    {
        "x": 0.9,
        "y": 0.9,
        "text": "Bottom right",
        "size": 20,
        "rotation": 0,
        "zoomable": true,
        "placement": "BOTTOM_RIGHT"
    }
]

*/