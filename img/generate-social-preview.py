#!/usr/bin/env python3
"""Generate the GitHub social preview image (1280x640) for CHILS."""

from PIL import Image, ImageDraw, ImageFont
import math
import os
import random

random.seed(42)

W, H = 1280, 640
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "social-preview.png")

canvas = Image.new("RGB", (W, H), (15, 23, 42))
draw = ImageDraw.Draw(canvas)

# Subtle radial gradient
cx, cy = W // 2, H // 2
for r in range(500, 0, -2):
    frac = r / 500
    c1 = int(15 + (25 - 15) * (1 - frac))
    c2 = int(23 + (35 - 23) * (1 - frac))
    c3 = int(42 + (62 - 42) * (1 - frac))
    draw.ellipse([cx - r * 1.6, cy - r, cx + r * 1.6, cy + r], fill=(c1, c2, c3))

bg = (15, 23, 42)

# Colors matching CHILS banner: teal (IS nodes) and red (non-IS nodes)
teal = (0, 178, 148)
red = (200, 60, 60)
dim_gray = (80, 90, 110)

# We'll show 3 solution graphs converging with arrows to a D-Core
# Each small graph has the same structure but different IS selections

# Base graph structure (small grid-like graph)
base_nodes = [
    (0, 0), (40, -30), (80, 0), (0, 50), (40, 25),
    (80, 50), (20, 90), (60, 90),
]

base_edges = [
    (0, 1), (1, 2), (0, 3), (1, 4), (2, 5),
    (3, 4), (4, 5), (3, 6), (4, 6), (4, 7), (5, 7), (6, 7),
    (0, 4), (1, 5), (2, 4),
]

# Different IS selections for each "solution"
solutions_is = [
    {0, 2, 6, 5},    # Solution 1
    {1, 3, 5, 7},    # Solution 2  -- different IS
    {0, 5, 6},        # Solution 3
]

# D-Core: the consensus (nodes that appear in most solutions are "strong")
dcore_strong = {5, 0}  # nodes in multiple IS solutions
dcore_weak = {1, 2, 3, 6, 7}

# Offsets for the 3 solution graphs
solution_offsets = [
    (100, 120),
    (100, 290),
    (100, 460),
]
dcore_offset = (380, 290)

node_r = 7
scale = 0.9

def draw_graph(offset_x, offset_y, is_set, label, label_y_offset=100):
    """Draw a small graph with highlighted IS nodes."""
    ox, oy = offset_x, offset_y
    # Edges
    for i, j in base_edges:
        x1 = int(ox + base_nodes[i][0] * scale)
        y1 = int(oy + base_nodes[i][1] * scale)
        x2 = int(ox + base_nodes[j][0] * scale)
        y2 = int(oy + base_nodes[j][1] * scale)
        draw.line([(x1, y1), (x2, y2)], fill=(50, 60, 80), width=1)

    # Nodes
    for idx, (nx, ny) in enumerate(base_nodes):
        x = int(ox + nx * scale)
        y = int(oy + ny * scale)
        if idx in is_set:
            # Glow
            for gr in range(16, node_r, -1):
                af = 1 - (gr - node_r) / (16 - node_r)
                gc = tuple(int(teal[k] * 0.12 * af + (1 - 0.12 * af) * bg[k]) for k in range(3))
                draw.ellipse([x - gr, y - gr, x + gr, y + gr], fill=gc)
            draw.ellipse([x - node_r, y - node_r, x + node_r, y + node_r], fill=teal,
                         outline=(0, 220, 180), width=1)
        else:
            draw.ellipse([x - node_r, y - node_r, x + node_r, y + node_r], fill=red,
                         outline=(240, 80, 80), width=1)


# Draw the 3 solutions
for idx, (ox, oy) in enumerate(solution_offsets):
    draw_graph(ox, oy, solutions_is[idx], f"Sol {idx+1}")

# Draw arrows from solutions to D-Core area
arrow_color = (100, 115, 140)
for ox, oy in solution_offsets:
    start_x = ox + 95
    start_y = oy + 40
    end_x = dcore_offset[0] - 15
    end_y = dcore_offset[1] + 40
    # Draw arrow line
    draw.line([(start_x, start_y), (end_x, end_y)], fill=arrow_color, width=2)
    # Arrowhead
    angle = math.atan2(end_y - start_y, end_x - start_x)
    ah_len = 10
    for da in [-0.4, 0.4]:
        ax = end_x - ah_len * math.cos(angle + da)
        ay = end_y - ah_len * math.sin(angle + da)
        draw.line([(end_x, end_y), (ax, ay)], fill=arrow_color, width=2)

# Draw D-Core graph
dcx, dcy = dcore_offset
for i, j in base_edges:
    x1 = int(dcx + base_nodes[i][0] * scale)
    y1 = int(dcy + base_nodes[i][1] * scale)
    x2 = int(dcx + base_nodes[j][0] * scale)
    y2 = int(dcy + base_nodes[j][1] * scale)
    # Edges to weak nodes are dashed/dim
    if i in dcore_weak or j in dcore_weak:
        draw.line([(x1, y1), (x2, y2)], fill=(40, 50, 65), width=1)
    else:
        draw.line([(x1, y1), (x2, y2)], fill=(60, 75, 100), width=2)

for idx, (nx, ny) in enumerate(base_nodes):
    x = int(dcx + nx * scale)
    y = int(dcy + ny * scale)
    if idx in dcore_strong:
        # Strong D-Core nodes: solid black with glow
        for gr in range(14, 6, -1):
            af = 1 - (gr - 6) / 8
            gc = tuple(int(220 * 0.1 * af + (1 - 0.1 * af) * bg[k]) for k in range(3))
            draw.ellipse([x - gr, y - gr, x + gr, y + gr], fill=gc)
        draw.ellipse([x - node_r, y - node_r, x + node_r, y + node_r], fill=(220, 225, 235),
                     outline=(255, 255, 255), width=2)
    else:
        # Weak/faded nodes
        draw.ellipse([x - node_r + 1, y - node_r + 1, x + node_r - 1, y + node_r - 1],
                     fill=(50, 60, 75), outline=(70, 80, 100), width=1)

# Labels under each graph
try:
    font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14)
except OSError:
    font_small = ImageFont.load_default()

for idx, (ox, oy) in enumerate(solution_offsets):
    draw.text((ox + 10, oy + 92), f"solution {idx+1}", fill=(100, 115, 140), font=font_small)
draw.text((dcx + 10, dcy + 92), "D-Core", fill=(160, 170, 190), font=font_small)

# Fonts for text side
try:
    font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 80)
    font_sub = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 28)
    font_tag = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
    font_legend = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 18)
except OSError:
    font_title = font_sub = font_tag = font_legend = ImageFont.load_default()

text_x = 580

# Title
draw.text((text_x, 160), "CHILS", fill=(240, 245, 255), font=font_title)

# Separator
draw.line([(text_x, 265), (text_x + 560, 265)], fill=(60, 80, 110), width=2)

# Subtitle
draw.text((text_x, 285), "Concurrent Hybrid Iterated", fill=(180, 195, 220), font=font_sub)
draw.text((text_x, 323), "Local Search", fill=(180, 195, 220), font=font_sub)

# Tagline
draw.text((text_x, 390), "Fast heuristic for the maximum weight", fill=(110, 130, 160), font=font_tag)
draw.text((text_x, 420), "independent set problem", fill=(110, 130, 160), font=font_tag)

# Badge: Part of KaMIS
badge_y = 490
badge_text = "Part of KaMIS"
badge_bbox = draw.textbbox((0, 0), badge_text, font=font_legend)
bw = badge_bbox[2] - badge_bbox[0]
draw.rounded_rectangle(
    [text_x, badge_y, text_x + bw + 24, badge_y + 32],
    radius=16, fill=(38, 50, 56)
)
draw.text((text_x + 12, badge_y + 5), badge_text, fill=(0, 200, 170), font=font_legend)

# Decorative dots
for i, dx in enumerate(range(0, 300, 40)):
    dot_x = text_x + dx
    dot_y = 560
    if i % 2 == 0:
        draw.ellipse([dot_x - 5, dot_y - 5, dot_x + 5, dot_y + 5], fill=teal)
    else:
        draw.ellipse([dot_x - 3, dot_y - 3, dot_x + 3, dot_y + 3], fill=(60, 80, 110))

canvas.save(OUT_PATH, "PNG", quality=95)
print(f"Saved {OUT_PATH}")
