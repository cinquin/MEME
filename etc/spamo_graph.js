"use strict";

var SpamoConst1 = {
  "size_title" : 18,
  "font_title" : "18px Helvetica",
  "size_subtitle" : 14,
  "font_subtitle" : "14px Helvetica",
  //font for the axis descriptions
  "size_axis" : 12,
  "font_axis" : "12px Helvetica",
  //font for the numbers on the axis
  "size_label" : 9,
  "font_label" : "9px Helvetica",
  "size_tic" : 2,
  "spacer_tic" : 2,
  "spacer_side" : 5,
  "spacer_top" : 5,
  "spacer_bottom" : 5,
  "spacer_title" : 0,
  "spacer_subtitle" : 5,
  "spacer_axis" : 0,
  "spacer_base" : 3,
  "spacer_mid" : 10,
  "spacer_border" : 5,
  "y_axis_label_start" : "Number of Occurrences (total = ",
  "y_axis_label_end" : ")",
  "x_axis_label" : "Distance from Primary to Secondary Motif (gap)",
  "same_strand_label" : "Same Strand",
  "oppo_strand_label" : "Opposite Strand",
  "upstream_label" : "Upstream",
  "downstream_label" : "Downstream",
};


/*
 * Measures the character heights in string
 * by painting them over the top of each other and finding the
 * bounding box. Returns the height of the total string
 * as well as calculating the vertical centering offset.
 */
function string_height(text, font, size) {
  var i, x, y, box;
  if (typeof text === "undefined") throw new Error("Text is not set!");
  if (typeof font === "undefined") throw new Error("Font is not set!");
  if (typeof size === "undefined") throw new Error("Font size is not set!");
  x = size / 2;
  y = size + (size / 2);
  if (!string_height.box) string_height.box = document.createElement("canvas");
  string_height.box.width = size * 2;
  string_height.box.height = size * 2;
  var ctx = string_height.box.getContext('2d');
  ctx.font = font;
  for (i = 0; i < text.length; ++i) {
    ctx.fillText(text.charAt(i), x, y);
  }
  box = bbox(ctx, size * 2, size * 2, true);

  return {height: box.height, vcenteroffset: (box.height / 2)}
}

/*
 * Calculates a tight bounding box on the contents of a canvas
 */
function bbox(ctx, cwidth, cheight, skip_width_calculation, skip_height_calculation) {
  if (typeof skip_width_calculation != 'boolean') skip_width_calculation = false;
  if (typeof skip_height_calculation != 'boolean') skip_height_calculation = false;
  if (skip_width_calculation && skip_height_calculation) throw new Error("bbox: can't skip calculation of both width and height");
  var data = ctx.getImageData(0, 0, cwidth, cheight).data;
  var r = 0, c = 0;// r: row, c: column
  var top_line = -1, bottom_line = -1, left_line = -1, right_line = -1;
  var box_width = 0, box_height = 0;
  if (!skip_height_calculation) {
    // Find the top-most line with a non-white pixel
    for (r = 0; r < cheight; r++) {
      for (c = 0; c < cwidth; c++) {
        if (data[r * cwidth * 4 + c * 4 + 3]) {
          top_line = r;
          break;
        }
      }
      if (top_line != -1) break;
    }
    
    //find the last line with a non-white pixel
    if (top_line != -1) {
      for (r = cheight-1; r >= top_line; r--) {
        for(c = 0; c < cwidth; c++) {
          if(data[r * cwidth * 4 + c * 4 + 3]) {
            bottom_line = r;
            break;
          }
        }
        if (bottom_line != -1) break;
      }
      box_height = bottom_line - top_line + 1;
    }
  }

  if (!skip_width_calculation) {
    // Find the left-most line with a non-white pixel
    for (c = 0; c < cwidth; c++) {
      for (r = 0; r < cheight; r++) {
        if (data[r * cwidth * 4 + c * 4 + 3]) {
          left_line = c;
          break;
        }
      }
      if (left_line != -1) break;
    }

    //find the right most line with a non-white pixel
    if (left_line != -1) {
      for (c = cwidth-1; c >= left_line; c--) {
        for(r = 0; r < cheight; r++) {
          if(data[r * cwidth * 4 + c * 4 + 3]) {
            right_line = c;
            break;
          }
        }
        if (right_line != -1) break;
      }
      box_width = right_line - left_line + 1;
    }
  }
  //return the bounds
  if (!skip_height_calculation && !skip_width_calculation) {
    return {top: top_line, bottom: bottom_line, left: left_line, right: right_line, width: box_width, height: box_height};
  } else if (!skip_height_calculation) {
    return {top: top_line, bottom: bottom_line, height: box_height};
  } else {
    return {left: left_line, right: right_line, width: box_width};
  }
}

var SpamoQuadGraphMetrics = function(spamo, target_width, target_height, ctx) {
  var i;

  this.width = target_width;
  this.height = target_height;

  this.revcomp = (spamo.bins.length == 4);

  //count axis numbers width
  this.count_width = 0;
  ctx.font = SpamoConst1.font_label;
  for (i = 0; i <= spamo.y_axis_max; i+= 5) {
    var width = ctx.measureText("" + i).width;
    if (width > this.count_width) {
      this.count_width = width;
    }
  }
  //spacing axis numbers width
  this.spacing_width = 0;
  ctx.font = SpamoConst1.font_label;
  for (i = -spamo.graph_margin; i <= spamo.graph_margin; ++i) {
    var width = ctx.measureText("" + i).width;
    if (width > this.spacing_width) {
      this.spacing_width = width;
    }
  }
  //calculate the height of a quadrant
  this.quad_height = this.height;
  if (spamo.graph_has_border) {
    this.quad_height -= 2 * SpamoConst1.spacer_border;
  }
  this.quad_height -= SpamoConst1.spacer_top;
  if (spamo.graph_has_title) {
    this.quad_height -= SpamoConst1.size_title;
  }
  if (spamo.graph_has_title && spamo.graph_has_subtitle) {
    this.quad_height -= SpamoConst1.spacer_title;
  }
  if (spamo.graph_has_subtitle) {
    this.quad_height -= SpamoConst1.size_subtitle;
  }
  if (spamo.graph_has_title || spamo.graph_has_subtitle) {
    this.quad_height -= SpamoConst1.spacer_subtitle;
  }
  if (spamo.graph_has_stream_label) {
    this.quad_height -= SpamoConst1.size_axis;
  }
  this.quad_height -= SpamoConst1.spacer_base * (this.revcomp ? 2 : 1);
  this.quad_height -= this.spacing_width;
  if (spamo.graph_has_x_axis_label) {
    this.quad_height -= SpamoConst1.spacer_axis + SpamoConst1.size_axis;
  }
  this.quad_height -= SpamoConst1.spacer_bottom;
  if (this.revcomp) this.quad_height /= 2;
  this.quad_height = Math.floor(this.quad_height);

  //calculate the max width of a quadrant
  this.quad_max_width = this.width;
  if (spamo.graph_has_border) {
    this.quad_max_width -= 2 * SpamoConst1.spacer_border;
  }
  this.quad_max_width -= 2 * SpamoConst1.spacer_side;
  if (spamo.graph_has_y_axis_label) {
    this.quad_max_width -= SpamoConst1.size_axis + SpamoConst1.spacer_axis; 
  }
  this.quad_max_width -= this.count_width + SpamoConst1.spacer_tic + SpamoConst1.spacer_mid;
  if (spamo.graph_has_strand_label) {
    this.quad_max_width -= SpamoConst1.size_axis + SpamoConst1.spacer_axis;
  }
  this.quad_max_width /= 2;
  this.quad_max_width = Math.floor(this.quad_max_width);

  //calculate the width of a graph bar the size of a single unit
  this.unit_width = this.quad_max_width / spamo.graph_margin;
  if (this.unit_width > 1) this.unit_width = Math.floor(this.unit_width);

  //calculate the width of the bar which represents the bin furtherest from the
  //primary which may be smaller than the other bars because the motif width
  //does not allow 
  this.partial_bar_width = ((spamo.graph_margin - spamo.graph_mlength + 1) % spamo.graph_binsize) * this.unit_width;

  //calculate the width of a full bar
  this.full_bar_width = this.unit_width * spamo.graph_binsize;

  //calculate the true width of a quadrant
  this.quad_width = this.unit_width * spamo.graph_margin;

  //add the excess onto the central spacer
  this.spacer_mid = SpamoConst1.spacer_mid + (this.quad_max_width - this.quad_width) * 2;

  //calculate the left quad
  this.left_quad = 0;
  if (spamo.graph_has_border) {
    this.left_quad += SpamoConst1.spacer_border;
  }
  this.left_quad += SpamoConst1.spacer_side;
  if (spamo.graph_has_y_axis_label) {
    this.left_quad += SpamoConst1.size_axis + SpamoConst1.spacer_axis;
  }
  this.left_quad += this.count_width + SpamoConst1.spacer_tic + SpamoConst1.size_tic;

  //calculate the right quad
  this.right_quad = this.left_quad + this.quad_width + this.spacer_mid;

  //calculate the bottom quad
  this.bottom_quad = 0;
  if (spamo.graph_has_border) {
    this.bottom_quad += SpamoConst1.spacer_border;
  }
  this.bottom_quad += SpamoConst1.spacer_bottom;
  if (spamo.graph_has_x_axis_label) {
    this.bottom_quad += SpamoConst1.size_axis + SpamoConst1.spacer_axis;
  }

  //calculate the top quad
  if (this.revcomp) {
    this.top_quad = this.bottom_quad + this.quad_height + (2 * SpamoConst1.spacer_base) + this.spacing_width;
  } else {
    this.top_quad = this.bottom_quad + this.spacing_width +  SpamoConst1.spacer_base;
  }
};

var SpamoQuadGraph = function (margin, bin_size, bin_max_count, secondary_motif_length, counts, highlights) {
  "use strict";
  var add = function (a, b) { return a + b; }
  //initialise graph state to default
  this.graph_title = "";
  this.graph_margin = margin;
  this.graph_mlength = secondary_motif_length
  this.graph_binsize = bin_size; 
  this.graph_binmax = bin_max_count;
  this.graph_sequences = counts.map(function(list) { return list.reduce(add); } ).reduce(add);
  this.graph_significant = 0.05;
  this.graph_has_border = false;
  this.graph_has_subtitle = false;
  this.graph_has_y_axis_label = true;
  this.graph_has_x_axis_label = true;
  this.graph_has_stream_label = false;
  this.graph_has_strand_label = false;

  //calculate bin count to read in the bin values
  var possible_count = this.graph_margin - this.graph_mlength + 1;
  var bin_count = Math.floor(possible_count / this.graph_binsize) + (possible_count % this.graph_binsize > 0 ? 1 : 0);
  var i;

  //get the bins
  this.bins = counts;
  this.highlights = highlights;

  //calculate variables
  this.graph_has_title = this.graph_title.length > 0;
  this.graph_subtitle = "margin: " + this.graph_margin + "   bin size: " + this.graph_binsize;
  this.y_axis_max = (Math.floor(this.graph_binmax / 5) + 1) * 5;
  this.full_bin_count = Math.floor((this.graph_margin - this.graph_mlength + 1) / this.graph_binsize);
};

SpamoQuadGraph.OPPOSITE_RIGHT = 3;
SpamoQuadGraph.OPPOSITE_LEFT = 2;
SpamoQuadGraph.SAME_RIGHT = 1;
SpamoQuadGraph.SAME_LEFT = 0;


SpamoQuadGraph.prototype.draw_border = function(metrics, ctx) {
  var bdr = SpamoConst1.spacer_border;
  var bdr2 = 2 * bdr;
  //add on 0.5 so the rectangle is in the center of a pixel and so is only 1 pixel wide
  ctx.strokeRect(bdr + 0.5, bdr + 0.5, metrics.width - bdr2, metrics.height - bdr2); 
};

SpamoQuadGraph.prototype.draw_title = function(metrics, ctx) {
  ctx.font = SpamoConst1.font_title;
  var title_width = ctx.measureText(this.graph_title).width;
  var x = (metrics.width / 2) - (title_width / 2);
  var y = metrics.height - SpamoConst1.spacer_top - SpamoConst1.size_title;
  if (this.graph_has_border) {
    y -= SpamoConst1.spacer_border;
  }
  ctx.fillText(this.graph_title, x, metrics.height - y);
};

SpamoQuadGraph.prototype.draw_subtitle = function(metrics, ctx) {
  ctx.font = SpamoConst1.font_subtitle;
  var subtitle_width = ctx.measureText(this.graph_subtitle).width;
  var x = (metrics.width / 2) - (subtitle_width / 2);
  var y = metrics.height - SpamoConst1.spacer_top - SpamoConst1.size_subtitle;
  if (this.graph_has_border) {
    y -= SpamoConst1.spacer_border;
  }
  if (this.graph_has_title) {
    y -= (SpamoConst1.size_title + SpamoConst1.spacer_title);
  }
  ctx.fillText(this.graph_subtitle, x, metrics.height - y);
};

SpamoQuadGraph.prototype.draw_y_axis_label = function(metrics, ctx) {
  ctx.font = SpamoConst1.font_axis;
  var txt = SpamoConst1.y_axis_label_start + this.graph_sequences + SpamoConst1.y_axis_label_end;
  var axislabel_width = ctx.measureText(txt).width;
  var x = SpamoConst1.spacer_side + SpamoConst1.size_axis;
  if (this.graph_has_border) {
    x += SpamoConst1.spacer_border;
  }
  var y = (metrics.height / 2) - (axislabel_width / 2);
  ctx.save();
  ctx.translate(x, metrics.height - y);
  ctx.rotate(-(Math.PI / 2));
  ctx.fillText(txt, 0, 0);
  ctx.restore();
};

SpamoQuadGraph.prototype.draw_x_axis_label = function(metrics, ctx) {
  ctx.font = SpamoConst1.font_axis;
  var axislabel_width = ctx.measureText(SpamoConst1.x_axis_label).width;
  var x = (metrics.width / 2) - (axislabel_width / 2);
  var y = SpamoConst1.spacer_bottom;
  if (this.graph_has_border) {
    y += SpamoConst1.spacer_border;
  }
  ctx.fillText(SpamoConst1.x_axis_label, x, metrics.height - y);
};

SpamoQuadGraph.prototype.draw_stream_labels = function(metrics, ctx) {
  ctx.font = SpamoConst1.font_axis;
  var upstream_width = ctx.measureText(SpamoConst1.upstream_label).width;
  var downstream_width = ctx.measureText(SpamoConst1.downstream_label).width;
  var y = metrics.top_quad + metrics.quad_height;
  var x = metrics.left_quad + (metrics.quad_width / 2) - (upstream_width / 2);
  ctx.fillText(SpamoConst1.upstream_label, x, metrics.height - y);
  x = metrics.right_quad + (metrics.quad_width / 2) - (downstream_width / 2);
  ctx.fillText(SpamoConst1.downstream_label, x, metrics.height - y);
};

SpamoQuadGraph.prototype.draw_strand_labels = function(metrics, ctx) {
  ctx.font = SpamoConst1.font_axis;
  var same_strand_width = ctx.measureText(SpamoConst1.same_strand_label).width;
  var x = metrics.right_quad + metrics.quad_width + SpamoConst1.spacer_axis + SpamoConst1.size_axis;
  var y = metrics.top_quad + (metrics.quad_height / 2) - (same_strand_width / 2);
  ctx.save();
  ctx.translate(x, metrics.height - y);
  ctx.rotate(-(Math.PI / 2));
  ctx.fillText(SpamoConst1.same_strand_label, 0, 0);
  ctx.restore();
  if (metrics.revcomp) {
    var oppo_strand_width = ctx.measureText(SpamoConst1.oppo_strand_label).width;
    y = metrics.bottom_quad + (metrics.quad_height / 2) - (oppo_strand_width / 2);
    ctx.save();
    ctx.translate(x, metrics.height - y);
    ctx.rotate(-(Math.PI / 2));
    ctx.fillText(SpamoConst1.oppo_strand_label, 0, 0);
    ctx.restore();
  }
};

SpamoQuadGraph.prototype.draw_y_axis_part = function(metrics, ctx, start_y, height) {
  ctx.font = SpamoConst1.font_label;
  ctx.beginPath();
  ctx.moveTo(metrics.left_quad - 0.5, metrics.height - start_y);
  ctx.lineTo(metrics.left_quad - 0.5, metrics.height - (start_y + height));
  ctx.stroke();
  //y_inc = -(height / this.y_axis_max * 5);
  var n_y_tics = 10
  var i_inc = Math.max(1, Math.floor(this.y_axis_max/n_y_tics));
  var y_inc = -(height / this.y_axis_max * i_inc);
  ctx.save();
  ctx.translate(metrics.left_quad, metrics.height - start_y);
  //for (i = 0; i <= this.y_axis_max; i += 5) {
  for (var i = 0; i <= this.y_axis_max; i += i_inc) {
    var lbl = "" + i;
    var lbl_width = ctx.measureText(lbl).width;
    var lbl_heightbox = string_height(lbl, SpamoConst1.font_label, SpamoConst1.size_label);
    var x = -(lbl_width + SpamoConst1.size_tic);
    var y = lbl_heightbox.vcenteroffset;
    //draw text
    ctx.fillText(lbl, x, y);

    //draw tic
    ctx.beginPath();
    ctx.moveTo(0, 0.5);
    ctx.lineTo(-SpamoConst1.size_tic, 0.5);
    ctx.stroke();

    ctx.translate(0, y_inc);
  }
  ctx.restore();
};

SpamoQuadGraph.prototype.draw_y_axis = function(metrics, ctx) {
  this.draw_y_axis_part(metrics, ctx, metrics.top_quad, metrics.quad_height);  
  if (metrics.revcomp) {
    this.draw_y_axis_part(metrics, ctx, metrics.bottom_quad + metrics.quad_height, -metrics.quad_height);  
  }
};

SpamoQuadGraph.prototype.draw_divider_part = function(metrics, ctx, x, y, len, dir) {
  var offset = 0;
  while (offset <= len) {
    ctx.fillRect(x, y + (offset * dir), 1, 2 * dir);
    offset += 3;
  }
};

SpamoQuadGraph.prototype.draw_divider = function(metrics, ctx) {
  this.draw_divider_part(metrics, ctx, metrics.left_quad + metrics.quad_width, 
      metrics.height - metrics.top_quad, metrics.quad_height, -1);

  this.draw_divider_part(metrics, ctx, metrics.right_quad - 1, 
      metrics.height - metrics.top_quad, metrics.quad_height, -1);

  if (metrics.revcomp) {
    this.draw_divider_part(metrics, ctx, metrics.left_quad + metrics.quad_width, 
        metrics.height - (metrics.bottom_quad + metrics.quad_height), metrics.quad_height, 1);

    this.draw_divider_part(metrics, ctx, metrics.right_quad - 1, 
        metrics.height - (metrics.bottom_quad + metrics.quad_height), metrics.quad_height, 1);
  }
};

SpamoQuadGraph.prototype.draw_motif_boundary = function(metrics, ctx) {
  if (this.graph_mlength == 1) return;
  var offset = metrics.unit_width * (this.graph_mlength - 1);

  this.draw_divider_part(metrics, ctx, metrics.left_quad + offset - 1,
      metrics.height - metrics.top_quad, metrics.quad_height, -1);

  this.draw_divider_part(metrics, ctx, metrics.right_quad + metrics.quad_width - offset, 
      metrics.height - metrics.top_quad, metrics.quad_height, -1);

  if (metrics.revcomp) {
    this.draw_divider_part(metrics, ctx, metrics.left_quad + offset - 1, 
        metrics.height - (metrics.bottom_quad + metrics.quad_height), metrics.quad_height, 1);

    this.draw_divider_part(metrics, ctx, metrics.right_quad + metrics.quad_width - offset, 
        metrics.height - (metrics.bottom_quad + metrics.quad_height), metrics.quad_height, 1);
  }
};

SpamoQuadGraph.prototype.draw_ender = function(metrics, ctx) {
  var left = metrics.right_quad + metrics.quad_width + 0.5;
  ctx.save();
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(left, metrics.height - metrics.top_quad);
  ctx.lineTo(left, metrics.height - (metrics.top_quad + metrics.quad_height));
  ctx.stroke();

  if (metrics.revcomp) {
    ctx.beginPath();
    ctx.moveTo(left, metrics.height - (metrics.bottom_quad + metrics.quad_height));
    ctx.lineTo(left, metrics.height - metrics.bottom_quad);
    ctx.stroke();
  }

  ctx.restore();
};

SpamoQuadGraph.prototype.draw_x_axis_num = function(metrics, ctx, num) {
  ctx.font = SpamoConst1.font_label;
  var txt = "" + num;
  var width = ctx.measureText(txt).width;
  var vcenter = string_height(txt, SpamoConst1.font_label, SpamoConst1.size_label).vcenteroffset;
  ctx.save();
  ctx.rotate(-(Math.PI / 2));
  ctx.fillText(txt, -width - 0.5, vcenter);
  ctx.restore();
};

SpamoQuadGraph.prototype.draw_x_axis = function(metrics, ctx) {
  var x = metrics.left_quad + metrics.quad_width;
  var y = metrics.top_quad - SpamoConst1.spacer_base;
  var skip = Math.ceil(SpamoConst1.size_label / metrics.unit_width);
  if (skip == 0) skip = 1;
  if (skip % 10 > 0) {
    skip -= (skip % 10) + 10;
  }
  ctx.save();
  ctx.translate(x - 0.5, metrics.height - y);
  for (var num = 0; num < this.graph_margin; num++) {
    if (num % skip == 0) {
      this.draw_x_axis_num(metrics, ctx, num);
    }
    ctx.translate(-metrics.unit_width, 0);
  }
  ctx.restore();
  x = metrics.right_quad;
  ctx.save();
  ctx.translate(x + 0.5, metrics.height - y);
  for (var num = 0; num < this.graph_margin; num++) {
    if (num % skip == 0) {
      this.draw_x_axis_num(metrics, ctx, num);
    }
    ctx.translate(metrics.unit_width, 0);
  }
  ctx.restore();

};

SpamoQuadGraph.prototype.draw_bar = function(metrics, ctx, bar_width, bin, highlight) {
  var bar_height;
  if (typeof bin !== "number" || bin == 0) return;
  bar_height = (bin / this.y_axis_max) * metrics.quad_height;
  ctx.fillStyle = (typeof highlight === "string" ? highlight : "rgb(205, 205, 205)");
  ctx.fillRect(0, 0, bar_width, -bar_height);
};

SpamoQuadGraph.prototype.draw_bars_part = function(metrics, ctx, bins, highlights) {
  var i;
  //skip unused
  if (this.graph_mlength > 1) ctx.translate(metrics.unit_width * (this.graph_mlength - 1), 0);
  i = Math.ceil((this.graph_margin - this.graph_mlength + 1) / this.graph_binsize) - 1;
  //draw partial bar
  if (metrics.partial_bar_width > 0) {
    this.draw_bar(metrics, ctx, metrics.partial_bar_width, bins[i], highlights[i]);
    ctx.translate(metrics.partial_bar_width, 0);
    i--;
  }
  //draw full bar
  for (; i >= 0; i--) {
    this.draw_bar(metrics, ctx, metrics.full_bar_width, bins[i], highlights[i]);
    ctx.translate(metrics.full_bar_width, 0);
  }
};

SpamoQuadGraph.prototype.draw_bars = function(metrics, ctx) {
  if (metrics.revcomp) {
    // bottom right
    ctx.save();
    ctx.translate(metrics.right_quad + metrics.quad_width, metrics.height - (metrics.bottom_quad + metrics.quad_height));
    ctx.scale(-1, -1);
    this.draw_bars_part(metrics, ctx, this.bins[SpamoQuadGraph.OPPOSITE_RIGHT], this.highlights[SpamoQuadGraph.OPPOSITE_RIGHT]);
    ctx.restore();
    // bottom left
    ctx.save();
    ctx.translate(metrics.left_quad, metrics.height - (metrics.bottom_quad + metrics.quad_height));
    ctx.scale(1, -1);
    this.draw_bars_part(metrics, ctx, this.bins[SpamoQuadGraph.OPPOSITE_LEFT], this.highlights[SpamoQuadGraph.OPPOSITE_LEFT]);
    ctx.restore();
  }
  // top right
  ctx.save();
  ctx.translate(metrics.right_quad + metrics.quad_width, metrics.height - metrics.top_quad);
  ctx.scale(-1, 1);
  this.draw_bars_part(metrics, ctx, this.bins[SpamoQuadGraph.SAME_RIGHT], this.highlights[SpamoQuadGraph.SAME_RIGHT]);
  ctx.restore();
  // top left
  ctx.save();
  ctx.translate(metrics.left_quad, metrics.height - metrics.top_quad);
  this.draw_bars_part(metrics, ctx, this.bins[SpamoQuadGraph.SAME_LEFT], this.highlights[SpamoQuadGraph.SAME_LEFT]);
  ctx.restore();
};

SpamoQuadGraph.prototype.draw_graph = function(ctx, w, h) {
  var metrics = new SpamoQuadGraphMetrics(this, w, h, ctx);
  
  if (this.graph_has_border) this.draw_border(metrics, ctx);
  if (this.graph_has_title) this.draw_title(metrics, ctx);
  if (this.graph_has_subtitle) this.draw_subtitle(metrics, ctx);
  if (this.graph_has_y_axis_label) this.draw_y_axis_label(metrics, ctx);
  if (this.graph_has_x_axis_label) this.draw_x_axis_label(metrics, ctx);
  if (this.graph_has_strand_label) this.draw_strand_labels(metrics, ctx);
  if (this.graph_has_stream_label) this.draw_stream_labels(metrics, ctx);
  this.draw_y_axis(metrics, ctx);
  this.draw_x_axis(metrics, ctx);
  this.draw_bars(metrics, ctx);
  this.draw_ender(metrics, ctx);
  this.draw_divider(metrics, ctx);
  this.draw_motif_boundary(metrics, ctx);
};

var SpamoConst2 = {
  "padding": 5,
  "xaxis_fontsize": 9,
  "xaxis_font": "9px Helvetica",
  "xaxis_fontcolour": "black",
  "xaxis_spacer": 2,
  "xlabel_fontsize": 12,
  "xlabel_font": "12px Helvetica",
  "xlabel_fontcolour": "black",
  "xlabel_spacer": 2,
  "yaxis_spacer": 1,
  "yaxis_tic": 3,
  "yaxis_ticspacer": 1,
  "yaxis_fontsize": 9,
  "yaxis_fontcolour": "black",
  "yaxis_font": "9px Helvetica",
  "yaxis_fontspacer": 4,
  "ylabel_fontsize": 12,
  "ylabel_fontcolour": "black",
  "ylabel_font": "12px Helvetica",
  "ylabel_spacer": 2,
  "bar_spacer": 1,
  "bar_colour": "rgb(205, 205, 205)"
};

var SpamoOrientGraphMetrics = function(graph, ctx, w, h) {
  "use strict";
  var i, txt_width;
  this.y_min = Math.floor(Math.max(graph.min_count - 1, 0) / 5) * 5;
  this.y_max = Math.ceil((graph.max_count + 1) / 5) * 5;
  // determine maximum width of y-axis text
  this.yaxis_width = 0;
  for (i = this.y_min; i <= this.y_max; i += 5) {
    txt_width = ctx.measureText("" + i).width;
    if (txt_width > this.yaxis_width) this.yaxis_width = txt_width;
  }

  this.width = w;
  this.height = h;

  this.graph_max_width = this.width - (SpamoConst2.padding * 2) -
    SpamoConst2.ylabel_fontsize - SpamoConst2.ylabel_spacer - this.yaxis_width -
    SpamoConst2.yaxis_ticspacer - SpamoConst2.yaxis_tic - SpamoConst2.yaxis_spacer;
  this.graph_max_height = this.height - (SpamoConst2.padding * 2) - 
    (SpamoConst2.yaxis_fontsize / 2) - SpamoConst2.xaxis_spacer - 
    SpamoConst2.xaxis_fontsize - SpamoConst2.xlabel_spacer -
    SpamoConst2.xlabel_fontsize;
  this.bar_width = Math.floor((this.graph_max_width / graph.range) - SpamoConst2.bar_spacer);
  this.bar_height = this.graph_max_height;
  this.graph_height = this.bar_height;
  this.graph_width = (this.bar_width + SpamoConst2.bar_spacer) * graph.range;
  // determine how many units the increment should be based on the font size and spacer size
  var y_units = this.y_max - this.y_min;
  var y_unit_height = this.bar_height / y_units;
  var y_font_units = (SpamoConst2.yaxis_fontsize + SpamoConst2.yaxis_fontspacer) / y_unit_height;
  this.y_inc = Math.ceil(y_font_units / 5) * 5; 
};

var SpamoOrientGraph = function(quadrant_counts, quadrant_list, start, range, bin_size, highlights) {
  "use strict";
  var i, j;
  var MAX_INT = 9007199254740992;
  // store the inputs
  this.hl = highlights;
  this.start = start; // starting bin
  this.range = range; // number of bins
  this.bin_size = bin_size; // size of bin
  // create the summed counts list
  this.counts = [];
  for (i = 0; i < range; i++) this.counts[i] = 0;
  for (i = 0; i < quadrant_list.length; i++) {
    for (j = 0; j < range; j++) {
      this.counts[j] += quadrant_counts[quadrant_list[i]][j + start];
    }
  }
  // calculate the minimum and maximum count
  this.max_count = 0;
  this.min_count = MAX_INT;
  for (i = 0; i < this.counts.length; i++) {
    if (this.counts[i] > this.max_count) this.max_count = this.counts[i];
    if (this.counts[i] < this.min_count) this.min_count = this.counts[i];
  }
};

SpamoOrientGraph.prototype.draw_bars = function(ctx, metrics) {
  var i, bar_height;
  ctx.save();
  for (i = 0; i < this.counts.length; i++) {
    bar_height = Math.round(((this.counts[i] - metrics.y_min) / (metrics.y_max - metrics.y_min)) * metrics.bar_height);
    ctx.fillStyle = (typeof this.hl[this.start + i] === "string" ? this.hl[this.start + i] : SpamoConst2.bar_colour);
    ctx.fillRect(0, 0, metrics.bar_width, -bar_height);
    ctx.translate(metrics.bar_width + SpamoConst2.bar_spacer, 0);
  }
  ctx.restore();
};

SpamoOrientGraph.prototype.draw_x_nums = function(ctx, metrics) {
  var i;
  ctx.save();
  ctx.font = SpamoConst2.xaxis_font;
  ctx.fillStyle = SpamoConst2.xaxis_fontcolour;
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  for (i = 0; i < this.counts.length; i++) {
    var col_text = "";
    if (this.bin_size == 1) {
      col_text = "" + (this.start + i);
    } else {
      var bin_start = (this.start + i) * this.bin_size;
      var bin_end = bin_start + this.bin_size - 1;
      col_text = "" + bin_start  + ":" + bin_end;
    }
    ctx.fillText(col_text, metrics.bar_width / 2, 0); 
    ctx.translate(metrics.bar_width + SpamoConst2.bar_spacer, 0);
  }
  ctx.restore();
};

SpamoOrientGraph.prototype.draw_y_nums = function(ctx, metrics) {
  var i, num, inc, y_pos;
  ctx.save();
  ctx.font = SpamoConst2.yaxis_font;
  ctx.fillStyle = SpamoConst2.yaxis_fontcolour;
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  ctx.lineWidth = 1;
  num = metrics.y_min;
  for (num = metrics.y_min; num <= metrics.y_max; num += metrics.y_inc) {
    y_pos = -( (num - metrics.y_min) / (metrics.y_max - metrics.y_min) ) * metrics.bar_height;
    ctx.fillText("" + num, 0, y_pos);
    ctx.beginPath();
    ctx.moveTo(SpamoConst2.yaxis_ticspacer, y_pos);
    ctx.lineTo(SpamoConst2.yaxis_ticspacer + SpamoConst2.yaxis_tic, y_pos);
    ctx.stroke();
  }
  ctx.fillRect(Math.round(SpamoConst2.yaxis_ticspacer + SpamoConst2.yaxis_tic) - 1, 0, 1, -metrics.bar_height);
  ctx.restore();
};

SpamoOrientGraph.prototype.draw_x_label = function(ctx, metrics) {
  ctx.save();
  ctx.font = SpamoConst2.xlabel_font;
  ctx.fillStyle = SpamoConst2.xlabel_fontcolour;
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.fillText("Distance from Primary to Secondary Motif (gap)", metrics.graph_width / 2, 0);
  ctx.restore();
};

SpamoOrientGraph.prototype.draw_y_label = function(ctx, metrics) {
  ctx.save();
  ctx.font = SpamoConst2.ylabel_font;
  ctx.fillStyle = SpamoConst2.ylabel_fontcolour;
  ctx.textAlign = "center";
  ctx.translate(0, metrics.graph_height / 2);
  ctx.rotate(-(Math.PI / 2));
  ctx.fillText("Number of Occurrences", 0, 0);
  ctx.restore();
};

SpamoOrientGraph.prototype.draw_graph = function(ctx, w, h) {
  "use strict";
  var metrics = new SpamoOrientGraphMetrics(this, ctx, w, h);
  ctx.save();
  ctx.translate(SpamoConst2.padding + SpamoConst2.ylabel_fontsize,
      SpamoConst2.padding + (SpamoConst2.yaxis_fontsize / 2));
  this.draw_y_label(ctx, metrics);
  ctx.translate(SpamoConst2.ylabel_spacer + metrics.yaxis_width,
       metrics.bar_height);
  this.draw_y_nums(ctx, metrics);
  ctx.translate(SpamoConst2.yaxis_ticspacer + SpamoConst2.yaxis_tic + SpamoConst2.yaxis_spacer, 0);
  this.draw_bars(ctx, metrics);
  ctx.save();
  ctx.translate(metrics.graph_width, 0);
  ctx.fillRect(0, 0, 1, -metrics.graph_height);
  ctx.restore();
  ctx.translate(0, SpamoConst2.xaxis_spacer);
  this.draw_x_nums(ctx, metrics);
  ctx.translate(0, SpamoConst2.xaxis_fontsize + SpamoConst2.xlabel_spacer);
  this.draw_x_label(ctx, metrics);
  ctx.restore();
};
