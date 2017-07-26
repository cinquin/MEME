
// PRIVATE GLOBAL (uhoh)
var _block_colour_lookup = {};

function block_colour(index) {
  function hsl2rgb(hue, saturation, lightness) {
    "use strict";
    function _hue(p, q, t) {
      "use strict";
      if (t < 0) t += 1;
      else if (t > 1) t -= 1;
      if (t < (1.0 / 6.0)) {
        return p + ((q - p) * 6.0 * t);
      } else if (t < 0.5) {
        return q;
      } else if (t < (2.0 / 3.0)) {
        return p + ((q - p) * ((2.0 / 3.0) - t) * 6.0);
      } else {
        return p;
      }
    }
    function _pad_hex(value) {
      var hex = Math.round(value * 255).toString(16);
      if (hex.length < 2) hex = "0" + hex;
      return hex;
    }
    var r, g, b, p, q;
    if (saturation == 0) {
      // achromatic (grayscale)
      r = lightness;
      g = lightness;
      b = lightness;
    } else {
      if (lightness < 0.5) {
        q = lightness * (1 + saturation);
      } else {
        q = lightness + saturation - (lightness * saturation);
      }
      p = (2 * lightness) - q;
      r = _hue(p, q, hue + (1.0 / 3.0));
      g = _hue(p, q, hue);
      b = _hue(p, q, hue - (1.0 / 3.0));
    }
    return "#" + _pad_hex(r) + _pad_hex(g) + _pad_hex(b);
  }
  if (typeof index !== "number" || index % 1 !== 0 || index < 0) return "#000000";
  // check for override
  if (_block_colour_lookup[index] == null) {
    var start = 0; //red
    var sat = 100;
    var light = 50;
    var divisions = 1 << Math.ceil(Math.log(index + 1) / Math.LN2);
    hue = start + (360 / divisions) * ((index - (divisions >> 1)) * 2 + 1);
    // colour input fields only support values in the form #RRGGBB
    _block_colour_lookup[index] = hsl2rgb(hue / 360, sat / 100, light / 100);
  }
  return _block_colour_lookup[index];
}

function set_block_colour(index, new_colour) {
  _block_colour_lookup[index] = new_colour;
  var blocks = document.querySelectorAll("div.block_motif[data-colour-index=\"" + index + "\"]");
  var i;
  for (i = 0; i < blocks.length; i++) {
    blocks[i].style.backgroundColor = new_colour;
  }
  var swatches = document.querySelectorAll("div.legend_swatch[data-colour-index=\"" + index + "\"]");
  var picker;
  for (i = 0; i < swatches.length; i++) {
    swatches[i].style.backgroundColor = new_colour;
    picker = swatches[i].querySelector("input[type=\"color\"]");
    if (picker != null) picker.value = new_colour;
  }
}

function make_block_legend_entry(motif_name, motif_colour_index) {
  if (typeof make_block_legend_entry.has_colour_picker !== "boolean") {
    // test if colour picker is supported, based off Modernizer
    // see http://stackoverflow.com/a/7787648/66387
    make_block_legend_entry.has_colour_picker = (function() {
      var doc_ele = document.documentElement;
      // We first check to see if the type we give it sticks..
      var input_ele = document.createElement('input');
      input_ele.setAttribute('type', 'color');
      var value_ok = input_ele.type !== 'text';
      if (value_ok) {
        // If the type does, we feed it a textual value, which shouldn't be valid.
        // If the value doesn't stick, we know there's input sanitization which infers a custom UI
        var smile = ':)';
        input_ele.value = smile;
        input_ele.style.cssText = 'position:absolute;visibility:hidden;';
        // chuck into DOM and force reflow for Opera bug in 11.00
        // github.com/Modernizr/Modernizr/issues#issue/159
        doc_ele.appendChild(input_ele);
        doc_ele.offsetWidth;
        value_ok = input_ele.value != smile;
        doc_ele.removeChild(input_ele);
      }
      return value_ok;
    })();
  }
  var entry = document.createElement("div");
  entry.className = "legend_entry";
  var swatch;
  swatch = document.createElement("div");
  swatch.className = "legend_swatch";
  swatch.setAttribute("data-colour-index", motif_colour_index);
  swatch.style.backgroundColor = block_colour(motif_colour_index);
  if (make_block_legend_entry.has_colour_picker) {
    var picker = document.createElement("input");
    picker.type = "color";
    picker.value = block_colour(motif_colour_index);
    picker.addEventListener("change", function(e) {
      set_block_colour(motif_colour_index, picker.value);
    }, false);
    swatch.addEventListener("click", function(e) {
      picker.click();
    }, false);
    swatch.appendChild(picker);
  }
  entry.appendChild(swatch);
  var name = document.createElement("div");
  name.className = "legend_text";
  name.appendChild(document.createTextNode(motif_name));
  entry.appendChild(name);
  return entry;
}

function make_block_ruler(max_len) {
  var container = document.createElement("div");
  container.className = "block_container";
  var step;
  if (max_len < 50) {
    step = 1;
  } else if (max_len < 100) {
    step = 2;
  } else if (max_len < 200) {
    step = 4;
  } else if (max_len < 500) {
    step = 10;
  } else if (max_len < 1000) {
    step = 20;
  } else if (max_len < 2000) {
    step = 40;
  } else if (max_len < 5000) {
    step = 100;
  } else if (max_len < 10000) {
    step = 200;
  } else if (max_len < 20000) {
    step = 400;
  } else {
    step = Math.floor(max_len / 20000) * 400;
  }
  var peroid;
  if (max_len < 10) {
    peroid = 1;
  } else if (max_len < 20) {
    peroid = 2;
  } else {
    peroid = 5;
  }
  var i, cycle, offset, tic, label;
  for (i = 0, cycle = 0; i < max_len; i += step, cycle = (cycle + 1) % peroid) {
    offset = "" + ((i / max_len) * 100) + "%";
    tic = document.createElement("div");
    tic.style.left = offset;
    tic.className = (cycle == 0 ? "tic_major" : "tic_minor");
    container.appendChild(tic);
    if (cycle == 0) {
      label = document.createElement("div");
      label.className = "tic_label";
      label.style.left = offset;
      label.appendChild(document.createTextNode(i));
      container.appendChild(label);
    }
  }
  return container;
}

function _calculate_block_needle_drag_pos(e, data) {
  var mouse;
  e = e || window.event;
  if (e.pageX || ev.pageY) {
    mouse = {"x": e.pageX, "y": e.pageY};
  } else {
    mouse = {
      x:e.clientX + document.body.scrollLeft - document.body.clientLeft, 
      y:e.clientY + document.body.scrollTop  - document.body.clientTop 
    };
  }
  var cont = data.container;
  var dragable_length = cont.clientWidth - 
    (cont.style.paddingLeft ? cont.style.paddingLeft : 0) -
    (cont.style.paddingRight ? cont.style.paddingRight : 0);
  //I believe that the offset parent is the body
  //otherwise I would need to make this recursive
  //maybe clientLeft would work, but the explanation of
  //it is hard to understand and it apparently doesn't work
  //in firefox 2.
  var diff = mouse.x - cont.offsetLeft;
  if (diff < 0) diff = 0;
  if (diff > dragable_length) diff = dragable_length;
  var pos = Math.round(diff / dragable_length * data.max);
  if (pos > data.len) pos = data.len;
  return pos;
}

function _update_block_needle_drag(e, data, done) {
  "use strict";
  var pos = _calculate_block_needle_drag_pos(e, data);
  // read the needle positions
  var left = parseInt(data.llabel.textContent, 10) - data.off - 1;
  var right = parseInt(data.rlabel.textContent, 10) - data.off;
  // validate needle positions
  if (left >= data.len) left = data.len - 1;
  if (left < 0) left = 0;
  if (right > data.len) right = data.len;
  if (right <= left) right = left + 1;
  // calculate the new needle positions
  if (data.moveboth) {
    var size = right - left;
    if (data.isleft) {
      if ((pos + size) > data.len) pos = data.len - size;
      left = pos;
      right = pos + size;
    } else {
      if ((pos - size) < 0) pos = size;
      left = pos - size;
      right = pos;
    }
  } else {
    if (data.isleft) {
      if (pos >= right) pos = right - 1;
      left = pos;
    } else {
      if (pos <= left) pos = left + 1;
      right = pos;
    }
  }
  // update the needle positions
  data.lneedle.style.left = "" + (left / data.max * 100) + "%";
  data.llabel.textContent = "" + (left + data.off + 1);
  data.rneedle.style.left = "" + (right / data.max * 100) + "%";
  data.rlabel.textContent = "" + (right + data.off);
  data.handler(left, right, done);
}

function _make_block_needle_drag_start_handler(isleft, data) {
  return function (e) {
    data.isleft = isleft;
    data.moveboth = !(e.shiftKey);
    document.addEventListener("mousemove", data.drag_during, false);
    document.addEventListener("mouseup", data.drag_end, false);
  };
}

function _make_block_needle_drag_end_handler(data) {
  return function (e) {
    document.removeEventListener("mousemove", data.drag_during, false);
    document.removeEventListener("mouseup", data.drag_end, false);
    _update_block_needle_drag(e, data, true);
  };
}

function _make_block_needle_drag_during_handler(data) {
  return function (e) {
    _update_block_needle_drag(e, data, false);
  };
}

// private function used by make_block_container
function _make_block_needle(isleft, value, data) {
  var vbar = document.createElement('div');
  vbar.className = "block_needle " + (isleft ? "left" : "right");
  vbar.style.left = "" + (value / data.max * 100)+ "%";
  var label = document.createElement('div');
  label.className = "block_handle " + (isleft ? "left" : "right");
  // The needles sit between the sequence positions, so the left one sits at the
  // start and the right at the end. This is why 1 is added to the displayed
  // value for a left handle as the user doesn't need to know about this detail
  label.textContent = "" + (isleft ? value + data.off + 1 : value + data.off);
  label.unselectable = "on"; // so IE and Opera don't select the text, others are done in css
  label.title = "Drag to move the displayed range. Hold shift and drag to change " + (isleft ? "lower" : "upper") + " bound of the range.";
  vbar.appendChild(label);
  if (isleft) {
    data.lneedle = vbar;
    data.llabel = label;
  } else {
    data.rneedle = vbar;
    data.rlabel = label;
  }
  label.addEventListener("mousedown", _make_block_needle_drag_start_handler(isleft, data), false);
  return vbar;
}

function make_block_container(is_stranded, has_both_strands, max_len, show_len, offset, range_handler) {
  offset = (offset != null ? offset : 0);
  // make the container for the block diagram
  var container = document.createElement("div");
  container.className = "block_container";
  container.setAttribute("data-max", max_len);
  container.setAttribute("data-off", offset);
  if (is_stranded) {
    var plus = document.createElement("div");
    plus.appendChild(document.createTextNode("+"));
    plus.className = "block_plus_sym";
    container.appendChild(plus);
    if (has_both_strands) {
      var minus = document.createElement("div");
      minus.appendChild(document.createTextNode("-"));
      minus.className = "block_minus_sym";
      container.appendChild(minus);
    }
  }
  var rule = document.createElement("div");
  rule.className = "block_rule";
  rule.style.width = ((show_len / max_len) * 100) + "%";
  container.appendChild(rule);
  if (range_handler != null) {
    var range_data = {
      "max": max_len,
      "len": show_len,
      "off": offset,
      "handler": range_handler,
      "container": container,
      "lneedle": null, "llabel": null,
      "rneedle": null, "rlabel": null,
      "isleft": false, "moveboth" : false
    };
    range_data.drag_during = _make_block_needle_drag_during_handler(range_data);
    range_data.drag_end = _make_block_needle_drag_end_handler(range_data);
    container.appendChild(_make_block_needle(false, 1, range_data)); // add right first so z-index works
    container.appendChild(_make_block_needle(true, 0, range_data));
  }
  return container;
}

function make_block_label(container, max_len, pos, length, message) {
  "use strict";
  var label = document.createElement("div");
  label.className = "block_label";
  label.style.left = (((pos + (length / 2)) / max_len) * 100) + "%";
  label.appendChild(document.createTextNode(message));
  container.appendChild(label);
}

function make_block(container, max_len,
    site_pos, site_len, site_pvalue, site_rc, site_colour_index, site_secondary) {
  "use strict";
  var block_height, block, block_region1, block_region2;
  var max_block_height = 12;
  var max_pvalue = 1e-10;
  // calculate the height of the block
  block_height = (site_pvalue < max_pvalue ? max_block_height : 
      (Math.log(site_pvalue) / Math.log(max_pvalue)) * max_block_height);
  if (block_height < 1) block_height = 1;
  // create a block to represent the motif
  block = document.createElement("div");
  block.className = "block_motif" + (site_secondary ? " scanned_site" : "") + (site_rc ? " bottom" : " top");
  block.style.left = ((site_pos / max_len) * 100) + "%";
  block.style.top = (!site_rc ? max_block_height - block_height : 
      max_block_height + 1) + "px";
  block.style.width = ((site_len / max_len) * 100) + "%";
  block.style.height = block_height + "px";
  block.style.backgroundColor = block_colour(site_colour_index);
  block.setAttribute("data-colour-index", site_colour_index);
  // add to container
  container.appendChild(block);
  var activator = function (e) {
    toggle_class(block, "active", true);
    var new_e = new e.constructor(e.type, e);
    block.dispatchEvent(new_e);
  };
  var deactivator = function (e) {
    toggle_class(block, "active", false);
    var new_e = new e.constructor(e.type, e);
    block.dispatchEvent(new_e);
  }
  // create a larger region to detect mouseover for the block
  block_region1 = document.createElement("div");
  block_region1.className = "block_region top" + 
    (site_secondary ? " scanned_site" : "") + (site_rc ? "" : " main");
  block_region1.style.left = block.style.left;
  block_region1.style.width = block.style.width;
  block_region1.addEventListener('mouseover', activator, false);
  block_region1.addEventListener('mouseout', deactivator, false);
  container.appendChild(block_region1);
  block_region2 = document.createElement("div");
  block_region2.className = "block_region bottom" + 
    (site_secondary ? " scanned_site" : "") + (site_rc ? " main" : "");
  block_region2.style.left = block.style.left;
  block_region2.style.width = block.style.width;
  block_region2.addEventListener('mouseover', activator, false);
  block_region2.addEventListener('mouseout', deactivator, false);
  container.appendChild(block_region2);
  return block;
}

function set_block_needle_positions(containingNode, start, end) {
  var container, lneedle, llabel, rneedle, rlabel, max, off, left, right;
  container = (/\bblock_container\b/.test(containingNode.className) ? containingNode : containingNode.querySelector(".block_container"));
  max = parseInt(container.getAttribute("data-max"), 10);
  off = parseInt(container.getAttribute("data-off"), 10);
  left = start - off;
  right = end - off;
  lneedle = containingNode.querySelector(".block_needle.left");
  llabel = lneedle.querySelector(".block_handle.left");
  rneedle = containingNode.querySelector(".block_needle.right");
  rlabel = rneedle.querySelector(".block_handle.right");
  // update the needle positions
  lneedle.style.left = "" + (left / max * 100) + "%";
  llabel.textContent = "" + (left + off + 1);
  rneedle.style.left = "" + (right / max * 100) + "%";
  rlabel.textContent = "" + (right + off);
}

function get_block_needle_positions(containingNode) {
  var container, llabel, rlabel, max, off, left, right;
  container = (/\bblock_container\b/.test(containingNode.className) ? containingNode : containingNode.querySelector(".block_container"));
  max = parseInt(container.getAttribute("data-max"), 10);
  off = parseInt(container.getAttribute("data-off"), 10);
  llabel = containingNode.querySelector(".block_needle.left > .block_handle.left");
  rlabel = containingNode.querySelector(".block_needle.right > .block_handle.right");
  left = parseInt(llabel.textContent, 10) - off - 1;
  right = parseInt(rlabel.textContent, 10) - off;
  return {"start": left + off, "end": right + off};
}
