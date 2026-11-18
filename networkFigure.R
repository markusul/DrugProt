# Static figures of the summary graph and the temporal graph of a query, in the style of Figure 3 of the paper.
# place_labels() and net_plot() follow R/pvalues/pValVis.R of https://github.com/markusul/DrugProt-Paper.
# static_network_figure() writes one of the two graphs as PNG.
library(ggplot2)

# a similar font is used if Arial is not installed
fig_family <- "Arial"
# node colors, separable under color blindness
fig_colors <- c(blue = "#2a78d6", orange = "#eb6834", green = "#1baf7a")
# edge colors for positive and negative dependencies in the temporal graph
sign_colors <- c(positive = "#4477aa", negative = "#cc4455")

# positions of the node labels
# each label is put at the candidate position around its node that covers no node and no other label
# and that crosses as little of the edges as possible, a label that is not next to its node gets a line to it
# nd and ed hold the nodes and edges in inches, the panel has the given width and height in inches
# the result is the same as in the paper, the cost of a position is only computed in two parts to save time
place_labels <- function(nd, ed, width, height, orders = 6) {
  lab <- which(nd$name != "")
  if (length(lab) == 0) {
    return(data.frame(name = character(0), lx = numeric(0), ly = numeric(0),
                      x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
                      hw = numeric(0), hh = numeric(0), node = integer(0)))
  }
  # half of the width and height of the text, with some space around it
  hw <- systemfonts::string_width(nd$name[lab], family = fig_family, size = 8) / 72 / 2 + 0.04
  hh <- 0.045 + 0.025
  # candidate positions, given by the direction from the node to the label in degrees and their distance
  cand <- expand.grid(angle = seq(0, 345, by = 15), d = c(0.06, 0.1, 0.16, 0.24, 0.34, 0.46, 0.6, 0.8))
  ux <- cos(cand$angle * pi / 180)
  uy <- sin(cand$angle * pi / 180)
  # only a label at the shortest distance is drawn without a line to its node
  has_line <- cand$d > 0.06

  # length of the part of a segment (x, y, xend, yend) that lies inside a box (xmin, xmax, ymin, ymax)
  # for all pairs of a box and a segment
  inside <- function(seg, box) {
    n <- nrow(box)
    m <- nrow(seg)
    out <- matrix(0, n, m)
    if (n == 0 || m == 0) return(out)
    # only a segment that comes near a box can lie inside it, the other pairs keep the length zero
    near <- which(outer(box[, 1], pmax(seg[, 1], seg[, 3]), "<=") & outer(box[, 2], pmin(seg[, 1], seg[, 3]), ">=") &
                    outer(box[, 3], pmax(seg[, 2], seg[, 4]), "<=") & outer(box[, 4], pmin(seg[, 2], seg[, 4]), ">="),
                  arr.ind = TRUE)
    if (nrow(near) == 0) return(out)
    b <- box[near[, 1], , drop = FALSE]
    s <- seg[near[, 2], , drop = FALSE]
    dx <- s[, 3] - s[, 1]
    dy <- s[, 4] - s[, 2]
    p <- cbind(-dx, dx, -dy, dy)
    q <- cbind(s[, 1] - b[, 1], b[, 2] - s[, 1], s[, 2] - b[, 3], b[, 4] - s[, 2])
    t0 <- rep(0, nrow(near))
    t1 <- rep(1, nrow(near))
    for (k in 1:4) {
      enters <- which(p[, k] < 0)
      leaves <- which(p[, k] > 0)
      t0[enters] <- pmax(t0[enters], q[enters, k] / p[enters, k])
      t1[leaves] <- pmin(t1[leaves], q[leaves, k] / p[leaves, k])
      t1[p[, k] == 0 & q[, k] < 0] <- -1
    }
    out[near] <- pmax(t1 - t0, 0) * sqrt(dx^2 + dy^2)
    out
  }

  # do two segments cross, for all pairs of a segment in a and a segment in b
  # the segments in a all start at the same point
  crosses <- function(a, b) {
    left <- function(px, py, qx, qy, rx, ry) (ry - py) * (qx - px) > (qy - py) * (rx - px)
    one <- function(v) matrix(v, nrow(a), nrow(b), byrow = TRUE)
    start <- left(a[1, 1], a[1, 2], b[, 1], b[, 2], b[, 3], b[, 4])
    end <- left(a[, 3], a[, 4], one(b[, 1]), one(b[, 2]), one(b[, 3]), one(b[, 4]))
    one(start) != end &
      left(a[, 1], a[, 2], a[, 3], a[, 4], one(b[, 1]), one(b[, 2])) !=
      left(a[, 1], a[, 2], a[, 3], a[, 4], one(b[, 3]), one(b[, 4]))
  }

  # a very large graph uses 1500 evenly chosen edges for the cost, the labels cannot avoid all edges there anyway
  edges <- as.matrix(ed)
  if (nrow(edges) > 1500) edges <- edges[round(seq(1, nrow(edges), length.out = 1500)), , drop = FALSE]
  # the candidates of each label and the part of their cost that does not depend on the other labels
  fixed <- lapply(seq_along(lab), function(i) {
    node <- c(nd$x[lab[i]], nd$y[lab[i]])
    # boxes of all candidates, the line of a candidate ends at the point of its box closest to the node
    r <- cand$d + pmin(hw[i] / abs(ux), hh / abs(uy))
    b <- cbind(node[1] + r * ux - hw[i], node[1] + r * ux + hw[i], node[2] + r * uy - hh, node[2] + r * uy + hh)
    l <- cbind(node[1], node[2], pmin(pmax(node[1], b[, 1]), b[, 2]), pmin(pmax(node[2], b[, 3]), b[, 4]))
    # short distances and labels without a line are preferred, and the label stays inside the panel
    cost <- 12 * cand$d + 30 * pmax(cand$d - 0.3, 0) + 6 * has_line +
      1e4 * (b[, 1] < 0.05 | b[, 2] > width - 0.05 | b[, 3] < 0.05 | b[, 4] > height - 0.05)
    # the label of a node without edges goes above the node if there is space
    far <- sqrt(pmin((edges[, 1] - node[1])^2 + (edges[, 2] - node[2])^2,
                     (edges[, 3] - node[1])^2 + (edges[, 4] - node[2])^2)) > 0.06
    if (all(far)) cost <- cost + 5 * (cand$angle != 90)
    # only the nodes within reach of a candidate can add to its cost
    reach <- 0.8 + 2 * (hw[i] + hh) + 0.21
    near <- which(abs(nd$x - node[1]) < reach & abs(nd$y - node[2]) < reach)
    nx <- nd$x[near]
    ny <- nd$y[near]
    other <- near != lab[i]
    # the label covers no node
    dist <- sqrt(pmax(outer(b[, 1], nx, "-"), 0, -outer(b[, 2], nx, "-"))^2 +
                   pmax(outer(b[, 3], ny, "-"), 0, -outer(b[, 4], ny, "-"))^2)
    cost <- cost + 300 * rowSums(dist < 0.035)
    # no other node is close to the label, so that it is clear which node it belongs to
    # this matters most for a label without a line
    cost <- cost + rowSums(dist[, other, drop = FALSE] < ifelse(has_line, 0.08, 0.2)) * ifelse(has_line, 40, 100)
    # the label crosses as little of the edges as possible
    cost <- cost + 60 * rowSums(inside(edges, b))
    # the line of the label passes no other node and no edge of another node
    vx <- l[, 3] - node[1]
    vy <- l[, 4] - node[2]
    t <- pmin(pmax((outer(vx, nx - node[1]) + outer(vy, ny - node[2])) / (vx^2 + vy^2), 0), 1)
    dist <- sqrt((outer(rep(1, nrow(cand)), nx - node[1]) - t * vx)^2 +
                   (outer(rep(1, nrow(cand)), ny - node[2]) - t * vy)^2)
    cost <- cost + 150 * has_line * rowSums(dist[, other, drop = FALSE] < 0.06)
    cost <- cost + 40 * has_line * rowSums(crosses(l, edges[far, , drop = FALSE]))
    list(b = b, l = l, cost = cost)
  })

  # places the labels one after the other in the given order, and then twice again given all other labels
  # returns the boxes of the labels, their lines to the nodes, and the total cost
  place <- function(first) {
    box <- matrix(NA_real_, length(lab), 4)
    line <- matrix(NA_real_, length(lab), 4)
    total <- numeric(length(lab))
    for (i in rep(first, 3)) {
      b <- fixed[[i]]$b
      l <- fixed[[i]]$l
      # the label does not overlap the other labels and their lines, and its line passes no other label
      others <- setdiff(which(!is.na(box[, 1])), i)
      ob <- box[others, , drop = FALSE]
      ol <- line[others, , drop = FALSE]
      ol <- ol[!is.na(ol[, 1]), , drop = FALSE]
      cost <- fixed[[i]]$cost +
        1e3 * rowSums(outer(b[, 1], ob[, 2], "<") & outer(b[, 2], ob[, 1], ">") &
                        outer(b[, 3], ob[, 4], "<") & outer(b[, 4], ob[, 3], ">")) +
        300 * rowSums(inside(ol, b) > 0) +
        300 * has_line * colSums(inside(l, ob) > 0)

      best <- which.min(cost)
      box[i, ] <- b[best, ]
      line[i, ] <- if (has_line[best]) l[best, ] else NA
      total[i] <- cost[best]
    }
    list(box = box, line = line, total = sum(total))
  }

  # the result depends on the order of the labels, the best of up to six orders is used
  # the orders start in a crowded region, in an empty region, or on one of the four sides
  crowded <- sapply(lab, function(v) sum(sqrt((nd$x - nd$x[v])^2 + (nd$y - nd$y[v])^2) < 0.6))
  res <- lapply(list(-crowded, crowded, nd$x[lab], -nd$x[lab], nd$y[lab], -nd$y[lab])[seq_len(orders)],
                function(v) place(order(v)))
  res <- res[[which.min(sapply(res, function(r) r$total))]]
  data.frame(name = nd$name[lab], lx = (res$box[, 1] + res$box[, 2]) / 2, ly = (res$box[, 3] + res$box[, 4]) / 2,
             x = res$line[, 1], y = res$line[, 2], xend = res$line[, 3], yend = res$line[, 4],
             hw = hw, hh = hh, node = lab)
}

# static version of a graph with a fixed layout, so that the same query always gives the same figure
# width and height are the size of the panel in inches, they are needed to place the labels
# edge_colour is one color or one color per row of Links
# Nodes$rank gives the priority of the labels, a label that does not fit is left out, the lowest priority first
net_plot <- function(Links, Nodes, colors, seed, width, height, edge_colour = "grey70") {
  g <- igraph::graph_from_data_frame(data.frame(from = Links$source + 1, to = Links$target + 1),
                                     directed = TRUE, vertices = data.frame(id = seq_len(nrow(Nodes))))
  g <- igraph::simplify(g)
  isolated <- igraph::degree(g) == 0
  # scales to the unit interval, nodes that share a coordinate go to the middle
  unit_scale <- function(v) if (diff(range(v)) > 0) (v - min(v)) / diff(range(v)) else rep(0.5, length(v))
  set.seed(seed)
  xy <- matrix(NA_real_, nrow(Nodes), 2)
  if (any(!isolated)) {
    xy[!isolated, ] <- igraph::layout_components(igraph::induced_subgraph(g, which(!isolated)),
                                                 layout = igraph::layout_with_fr)
    xy[!isolated, ] <- apply(xy[!isolated, , drop = FALSE], 2, unit_scale)
  }
  # nodes without an edge are placed in a row below the graph
  if (any(isolated)) {
    xy[isolated, ] <- cbind(if (sum(isolated) == 1) 0.5 else seq(0.08, 0.92, length.out = sum(isolated)),
                            if (all(isolated)) 0.5 else -0.14)
  }
  # the coordinates are in inches from here on, with a margin around the graph
  margin <- 0.25
  xy[, 1] <- margin + unit_scale(xy[, 1]) * (width - 2 * margin)
  xy[, 2] <- margin + unit_scale(xy[, 2]) * (height - 2 * margin)
  nd <- data.frame(x = xy[, 1], y = xy[, 2], name = Nodes$name, group = factor(Nodes$group, levels = names(colors)))
  # edges end shortly before the node they point to, self-dependencies are not drawn in the summary graph
  ed <- data.frame(x = nd$x[Links$source + 1], y = nd$y[Links$source + 1],
                   xend = nd$x[Links$target + 1], yend = nd$y[Links$target + 1],
                   colour = rep_len(edge_colour, nrow(Links)))
  ed <- ed[ed$x != ed$xend | ed$y != ed$yend, ]
  len <- sqrt((ed$xend - ed$x)^2 + (ed$yend - ed$y)^2)
  gap <- pmin(0.04 / len, 0.4)
  ed$xend <- ed$xend - gap * (ed$xend - ed$x)
  ed$yend <- ed$yend - gap * (ed$yend - ed$y)
  # a large graph gets fewer attempts for the labels, so that the figure is ready within seconds
  work <- sum(nd$name != "") * (nrow(nd) + nrow(ed))
  labels <- place_labels(nd, ed[, 1:4], width, height, orders = if (work < 5e3) 6 else if (work < 2e4) 2 else 1)
  # a dense graph has no room for all labels, those that overlap another label or leave the panel are left out
  labels$rank <- if (is.null(Nodes$rank)) seq_len(nrow(labels)) else Nodes$rank[labels$node]
  labels <- labels[labels$lx - labels$hw > 0 & labels$lx + labels$hw < width &
                     labels$ly - labels$hh > 0 & labels$ly + labels$hh < height, ]
  repeat {
    overlap <- outer(labels$lx - labels$hw, labels$lx + labels$hw, "<") & outer(labels$lx + labels$hw, labels$lx - labels$hw, ">") &
      outer(labels$ly - labels$hh, labels$ly + labels$hh, "<") & outer(labels$ly + labels$hh, labels$ly - labels$hh, ">")
    diag(overlap) <- FALSE
    if (!any(overlap)) break
    crowded <- which(rowSums(overlap) > 0)
    labels <- labels[-crowded[which.max(labels$rank[crowded])], ]
  }
  # white copies of the labels, shifted in all directions, separate the text from the edges
  halo <- do.call(rbind, lapply(seq(0, 2 * pi, length.out = 17)[-1], function(a) {
    transform(labels, lx = lx + 0.012 * cos(a), ly = ly + 0.012 * sin(a))
  }))

  ggplot(nd, aes(x, y)) +
    geom_segment(data = ed, aes(xend = xend, yend = yend), colour = ed$colour, linewidth = 0.3,
                 arrow = arrow(length = unit(0.04, "in"), type = "closed")) +
    geom_segment(data = labels[!is.na(labels$x), ], aes(xend = xend, yend = yend), linewidth = 0.25) +
    geom_text(data = halo, aes(lx, ly, label = name), family = fig_family, size = 8 / .pt, colour = "white") +
    geom_point(aes(colour = group), size = 1.3) +
    geom_text(data = labels, aes(lx, ly, label = name), family = fig_family, size = 8 / .pt) +
    scale_colour_manual(NULL, values = colors, drop = FALSE) +
    coord_cartesian(xlim = c(0, width), ylim = c(0, height), expand = FALSE) +
    theme_void(base_size = 10, base_family = fig_family) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8),
          plot.caption = element_text(size = 8, hjust = 0.5))
}

# priority of the proteins for a label, 1 is the highest and NA means no label
# selected proteins come first, then the proteins with the most edges
# a graph with at most max_labels nodes shows all names
# a larger graph shows the selected proteins and the proteins connected to more than one of them
# Links and Nodes are those of the summary graph
label_rank <- function(Links, Nodes, max_labels = 40) {
  selected <- Nodes$group == "Selected"
  pairs <- unique(rbind(as.matrix(Links[, c("source", "target")]), as.matrix(Links[, c("target", "source")])) + 1)
  pairs <- pairs[pairs[, 1] != pairs[, 2], , drop = FALSE]
  n_selected <- tabulate(pairs[selected[pairs[, 2]], 1], nrow(Nodes))
  degree <- tabulate(pairs[, 1], nrow(Nodes))
  keep <- if (nrow(Nodes) <= max_labels) rep(TRUE, nrow(Nodes)) else selected | n_selected > 1
  rank <- rep(NA_integer_, nrow(Nodes))
  rank[order(!keep, !selected, -degree)[seq_len(sum(keep))]] <- seq_len(sum(keep))
  rank
}

# writes the summary graph or the temporal graph of a query as PNG
# the arguments are the links and nodes that the application passes to forceNetwork()
# both graphs are needed for either figure, because the labels of the temporal graph follow the summary graph
# direction holds the edge colors of the temporal graph in the application, "blue" for a positive
# and "red" for a negative dependency
# size is the size of the figure in inches and dpi its resolution for a graph with up to 150 nodes
# a larger graph gets a larger figure with the same text size, up to twice the width and height,
# and a lower resolution, so that the number of pixels stays the same
static_network_figure <- function(file, graph = c("summary", "temporal"), Links_sum, Nodes_sum, Links_temp, Nodes_temp,
                                  direction = NULL, size = c(6, 5), dpi = 600) {
  graph <- match.arg(graph)
  # each graph tries at most max_labels names, more do not fit
  # a node of the temporal graph has the priority of its protein in the summary graph
  max_labels <- 40
  Nodes_sum$rank <- label_rank(Links_sum, Nodes_sum, max_labels)
  Nodes_temp$rank <- Nodes_sum$rank[match(sub("_[0-9]+$", "", Nodes_temp$name), Nodes_sum$name)]
  if (nrow(Nodes_temp) <= max_labels) Nodes_temp$rank[is.na(Nodes_temp$rank)] <- nrow(Nodes_sum) + 1
  Nodes_temp$rank[rank(Nodes_temp$rank, ties.method = "first", na.last = "keep") > max_labels] <- NA
  Nodes_sum$rank[!is.na(Nodes_sum$rank) & Nodes_sum$rank > max_labels] <- NA
  Nodes_sum$name[is.na(Nodes_sum$rank)] <- ""
  Nodes_temp$name[is.na(Nodes_temp$rank)] <- ""

  if (graph == "summary") {
    Links <- Links_sum
    Nodes <- Nodes_sum
    colors <- c(Selected = fig_colors[["orange"]], Connected = fig_colors[["blue"]])
    edge_colour <- "grey70"
    note <- NULL
  } else {
    Links <- Links_temp
    Nodes <- Nodes_temp
    times <- unique(Nodes_temp$group)
    colors <- setNames(rep_len(unname(fig_colors), length(times)), times)
    # edge colors by the sign of the dependency
    edge_colour <- "grey70"
    note <- NULL
    if (length(direction) == nrow(Links_temp) && any(direction %in% c("blue", "red"))) {
      edge_colour <- ifelse(direction == "blue", sign_colors[["positive"]],
                            ifelse(direction == "red", sign_colors[["negative"]], "grey70"))
      note <- "Blue edges are positive and red edges are negative dependencies"
    }
  }

  grow <- min(2, max(1, sqrt(nrow(Nodes) / 150)))
  size <- size * grow
  # the legend and the margin below it take 0.38 in of the height and the note 0.17 in
  below <- 0.38 + 0.17 * !is.null(note)
  fig <- net_plot(Links, Nodes, colors, seed = 1, width = size[1], height = size[2] - below, edge_colour = edge_colour) +
    labs(caption = note) +
    theme(plot.margin = margin(b = 4))
  ragg::agg_png(file, width = size[1], height = size[2], units = "in", res = dpi / grow, background = "white")
  on.exit(dev.off())
  print(fig)
  invisible(file)
}
