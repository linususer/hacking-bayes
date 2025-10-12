import math
from collections import defaultdict
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np

def get_idx(tree_depth, k):
    return int(((tree_depth + 1) * tree_depth) / 2 + k)

def get_parents(tree_depth, k):
    parent1 = get_idx(tree_depth - 1, k) if (k <= tree_depth - 1 and tree_depth > 0) else None
    parent2 = get_idx(tree_depth - 1, k - 1) if (k - 1 >= 0 and tree_depth > 0) else None
    return (parent1, parent2)

def add_node(tree, idx, parents, event_seq, p1, p2):
    heads, tails = event_seq
    count = math.comb(heads + tails, heads)
    p1_prob: np.float128 = (p1 ** heads) * ((1 - p1) ** tails)
    p2_prob: np.float128 = (p2 ** heads) * ((1 - p2) ** tails)
    BF = p1_prob / p2_prob if p2_prob != 0 else None
    tree[idx] = {
        'name': f"H{heads}T{tails}",
        'parent': parents,
        'event_value': event_seq,
        'count': count,
        'p1': p1_prob,
        'p2': p2_prob,
        'BF': BF,
        'stopped': False,
        'decision': "indecisive" # either "indecisive", "p1", or "p2"
    }

def create_tree(coinflips, p1, p2):
    tree = dict()
    tree_depth = 0
    idx = get_idx(tree_depth, 0)
    add_node(tree, idx, (None, None), (0,0), p1, p2)
    tree_depth = 1
    for i in range(1, coinflips + 1):
        for k in range(tree_depth + 1):
            idx = get_idx(tree_depth, k)
            parents = get_parents(tree_depth, k)
            event_seq = (tree_depth - k, k)
            add_node(tree, idx, parents, event_seq, p1, p2)
        tree_depth += 1
    return tree

def print_tree(tree):
    stopped_nodes = 0
    for idx, node in sorted(tree.items()):
        if node['stopped']:
            stopped_nodes += 1
        print(f"Node {idx}: {node['name']}, Parents: {node['parent']}, "
              f"Event Value: {node['event_value']}, Count: {node['count']}, "
              f"P1: {node['p1']:.4f}, P2: {node['p2']:.4f}, BF: {node['BF']}, decision: {node['decision']}, Stopped: {node['stopped']}")
    print(f"Total nodes: {len(tree)}, Stopped nodes: {stopped_nodes}")

def apply_fixed_sample_size_test(tree, coinflips, bf_crit, bf_crit2=None):
    if bf_crit2 is None:
        bf_crit2 = 1 / bf_crit
    tree_depth = coinflips
    for k in range(tree_depth + 1):
        idx = get_idx(tree_depth, k)
        node = tree[idx]
        bf = node['BF']
        if bf is not None:
            node['stopped'] = True  # "stopped" if BF exceeds threshold
            if bf > bf_crit:
                node['decision'] = "p1"
            elif bf < bf_crit2:
                node['decision'] = "p2"
            else:
                node['decision'] = "indecisive"
    return tree


def apply_optional_stopping(tree, coinflips, bf_crit, bf_crit2=None):
    if bf_crit2 is None:
        bf_crit2 = 1 / bf_crit
    for node in tree.values():
        BF = node.get('BF')
        if BF is not None and (BF > bf_crit or BF < (bf_crit2) or coinflips == sum(node['event_value'])):
            parents = node['parent']
            p1_idx, p2_idx = parents
            if len(parents) == 1 or len(parents) == 2:
                # if parents are not stopped, we count
                p1_stopped = tree.get(parents[0], {}).get('stopped', False) if parents[0] is not None else False
                p2_stopped = tree.get(parents[1], {}).get('stopped', False) if parents[1] is not None else False
                # all parents are stopped
                if ((p1_stopped and p2_stopped) or (p1_stopped and parents[1] is None) or
                    (p2_stopped and parents[0] is None)):
                    node['count'] = 0
                elif (p1_stopped or p2_stopped and parents[1] is not None and parents[0] is not None):
                    if (p1_stopped):
                        node['count'] = tree[p2_idx]['count']
                    elif (p2_stopped):
                        node['count'] = tree[p1_idx]['count']
                else:
                    count1 = tree[parents[0]]['count'] if parents[0] is not None else 0
                    count2 = tree[parents[1]]['count'] if parents[1] is not None else 0
                    node['count'] = count1 + count2
            node['stopped'] = True
            if BF > bf_crit:
                node['decision'] = "p1"
            elif BF < (bf_crit2):
                node['decision'] = "p2"
            else:
                node['decision'] = "indecisive"
        else:
            parents = node['parent']
            if parents is not None and len(parents) == 2:
                p1_idx, p2_idx = parents
                if p1_idx is not None and p2_idx is not None:
                    p1_stopped = tree.get(p1_idx, {}).get('stopped', False)
                    p2_stopped = tree.get(p2_idx, {}).get('stopped', False)
                    if p1_stopped and p2_stopped:
                        node['stopped'] = True
                        node['count'] = 0
                    elif p1_stopped or p2_stopped:
                        if p1_stopped:
                                node['count'] = tree[p2_idx]['count']
                        elif p2_stopped:
                                node['count'] = tree[p1_idx]['count']
                    else:
                        node['p1'] /= node['count']
                        node['p2'] /= node['count']
                        node['count'] = tree[p1_idx]['count'] + tree[p2_idx]['count']
                        node['p1'] = node['count'] * node['p1']
                        node['p2'] = node['count'] * node['p2']
                        node['BF'] = node['p1'] / node['p2'] if node['p2'] != 0 else None
            
            else:
                node['decision'] = "indecisive"
                node['stopped'] = False
    return tree

def draw_decision_box(fig, tree, pos, nodes, prob="p1", layer=3, padding=0.4, dash="dot", label_p1="R", label_p2="A", label_ind="I"):
    """
    Draw a decision box around nodes at a specific layer with a specific decision.
    """
    layer_nodes = []
    for n in nodes:
        d = 0
        while ((d+1)*d)//2 <= n:
            d += 1
        d -= 1
        if d == layer and tree[n]['decision'] == prob:
            layer_nodes.append(n)

    if layer_nodes:
        x_coords = [pos[n][0] for n in layer_nodes]
        y_coords = [pos[n][1] for n in layer_nodes]

        x0, x1 = min(x_coords) - padding, max(x_coords) + padding
        y0, y1 = min(y_coords) - padding, max(y_coords) + padding + 0.25

        fig.add_shape(
            type="rect",
            x0=x0,
            x1=x1,
            y0=y0,
            y1=y1,
            line=dict(color="black", width=3, dash=dash),
            fillcolor="rgba(0,0,0,0)"
        )

        if prob == "p1":
            label = label_p1
        elif prob == "p2":
            label = label_p2
        else:
            label = label_ind

        fig.add_annotation(
            x=x1 + 0.15,
            y=y1 - 0.25,
            text=label,
            showarrow=False,
            font=dict(size=30, color='black')
        )

    return fig

import numpy as np

def add_path_edges(fig, tree, pos, path, color='#CD1076', width=4, shorten=0.19):
    """
    Given a path in the tree, color the edges of this path.
    """
    edge_x, edge_y = [], []
    for i in range(len(path) - 1):
        n1, n2 = path[i], path[i+1]
        if n1 in pos and n2 in pos:
            x0, y0 = pos[n1]
            x1, y1 = pos[n2]
            # Shorten the line at both ends
            dx, dy = x1 - x0, y1 - y0
            dist = np.sqrt(dx**2 + dy**2)
            if dist > 0:
                factor = shorten / dist
                x0_new = x0 + dx * factor
                y0_new = y0 + dy * factor
                x1_new = x1 - dx * factor
                y1_new = y1 - dy * factor
            else:
                x0_new, y0_new, x1_new, y1_new = x0, y0, x1, y1

            edge_x.extend([x0_new, x1_new, None])
            edge_y.extend([y0_new, y1_new, None])

    fig.add_trace(go.Scatter(
        x=edge_x,
        y=edge_y,
        mode='lines',
        line=dict(color=color, width=width),
        hoverinfo='none',
        name='Path'
    ))
    return fig




def plot_tree_plotly(tree):
    nodes = list(tree.keys())
    labels = [tree[n]['name'] for n in nodes]
    counts = [tree[n]['count'] for n in nodes]
    
    pos = {}
    max_depth = 0
    for idx in nodes:
        d = 0
        while ((d+1)*d)//2 <= idx:
            d += 1
        d -= 1
        k = idx - ((d+1)*d)//2
        pos[idx] = (k - d/2, -d)
        max_depth = max(max_depth, d)
    
    # Build edge coordinates
    edge_x = []
    edge_y = []
    for idx, node in tree.items():
        x0, y0 = pos[idx]
        for parent in node['parent']:
            if parent is not None and parent in pos:
                x1, y1 = pos[parent]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])
    
    # Build nodes
    node_x = [pos[n][0] for n in nodes]
    node_y = [pos[n][1] for n in nodes]
    # text for counts
    text_x = []
    text_y = []
    text_labels = []
    node_color = []
    x_offset = 0.25
    y_offset = 0.02
    for n in nodes:
        x, y = pos[n]
        if tree[n]['count'] == 0:
            node_color.append('black')
        elif tree[n]['stopped']:
            node_color.append('orange')
            # if tree[n]['decision'] == "p1":
            #     # Add independent text (right side of node)
            #     text_x.append(x + x_offset)      # Offset right
            #     text_y.append(y - y_offset)      # Slight Y adjustment
            #     text_labels.append(f"-")
            # elif tree[n]['decision'] == "p2":
            #     # Add independent text (right side of node)
            #     text_x.append(x + x_offset)      # Offset right
            #     text_y.append(y - y_offset)      # Slight Y adjustment
            #     text_labels.append(f"+")
            # else:
            #     text_x.append(x + x_offset)      # Offset right
            #     text_y.append(y - y_offset)      # Slight Y adjustment
            #     text_labels.append(f"?")
        else:
            node_color.append('#59B3E6')


    fig = go.Figure()

    fig.add_trace(go.Scatter(x=edge_x, y=edge_y,
                             mode='lines',
                             line=dict(color='black', width=2),
                             hoverinfo='none'
                            ))
    
    fig.add_trace(go.Scatter(x=node_x, y=node_y,
                             mode='markers+text',
                             marker=dict(color=node_color, size=30, line_width=1, line_color='black'),
                             text=labels,
                             textposition="top center",
                             hoverinfo='text',
                             hovertext=[f"Node {n}: {tree[n]['name']}<br>BF: {tree[n]['BF']}<br>Count: {tree[n]['count']}" for n in nodes]
                            ))
    
    fig.add_trace(go.Scatter(x=text_x, y=text_y,
                                mode='text',
                                text=text_labels,
                                textfont=dict(size=20, color='black'),
                                hoverinfo='none'
                                ))
    
    # comment out for no decision boxes
    fig = draw_decision_box(fig, tree, pos, nodes, prob="p1", layer=5, padding=0.3, dash = "solid")
    fig = draw_decision_box(fig, tree, pos, nodes, prob="p2", layer=5, padding=0.3, dash = "solid")
    fig = draw_decision_box(fig, tree, pos, nodes, prob="indecisive", layer=5, padding=0.3, dash = "solid")

    fig = draw_decision_box(fig, tree, pos, nodes, prob="p1", layer=2, padding=0.3, dash = "dash", label_p1="R'", label_p2="A'", label_ind="I'")
    fig = draw_decision_box(fig, tree, pos, nodes, prob="p2", layer=2, padding=0.3, dash = "dash", label_p1="R'", label_p2="A'", label_ind="I'")
    fig = draw_decision_box(fig, tree, pos, nodes, prob="indecisive", layer=2, padding=0.3, dash = "dash", label_p1="R'", label_p2="A'", label_ind="I'")

    my_path = [0, 1, 3, 7, 12, 18]
    # comment out if no path should be colored
    fig = add_path_edges(fig, tree, pos, my_path, width=4)
    fig.update_layout(
        title="",
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='white',
        height=600,
        width=900,
        margin=dict(l=20, r=20, t=40, b=20)
    )
    for i, count in enumerate(counts):
        if count > 0:
            fig.add_annotation(
                x=node_x[i], y=node_y[i],
                text=str(count),
                showarrow=False,
                font=dict(size=20, color='black'),
                align='center',
            )
    fig.show(method="external")
    return fig

# print table of stopped nodes
def print_table(tree):
    nodes = [n for n in tree if tree[n]['count'] > 0 and tree[n]['stopped']]
    print("Node\tName\tCount\tP1\tP2\tBF\tDecision")
    # Sort nodes by event value (heads, tails) and then by count
    tree = {n: tree[n] for n in nodes}
    tree = dict(sorted(tree.items(), key=lambda item: (item[1]['BF'], item[1]['event_value'][0], item[1]['event_value'][1], item[1]['count']), reverse=True))
    leafs = []
    for idx, node in tree.items():
        parents = ', '.join(str(p) for p in node['parent'] if p is not None)
        if node['count'] != 0 and node['stopped']:
            leafs.append(node)
            print(f"{idx}\t{node['name']}\t{node['count']}\t"
                f"{node['p1']:.8f}\t{node['p2']:.8f}\t{round(node['BF'],3)}\t{node['decision']}")
    avg_tree_depth_p1 = sum([(node['event_value'][0] + node['event_value'][1]) * node['count'] * node['p1'] for node in tree.values()])
    avg_tree_depth_p2 = sum([(node['event_value'][0] + node['event_value'][1]) * node['count'] * node['p2'] for node in tree.values()])
    max_depth = max(node['event_value'][0] + node['event_value'][1] for node in tree.values())
    print(f"\nAverage tree depth for P(x | p1): {avg_tree_depth_p1:.4f} from max depth {max_depth}")
    print(f"\nAverage tree depth for P(x | p2): {avg_tree_depth_p2:.4f} from max depth {max_depth}")
    alpha_error = sum(node['p1'] * node['count'] for node in leafs if node['decision'] == "p2") + sum(node['p1'] * node['count'] for node in leafs if node['decision'] == "indecisive")
    beta_error = sum(node['p2'] * node['count'] for node in leafs if node['decision'] == "p1") + sum(node['p2'] * node['count'] for node in leafs if node['decision'] == "indecisive")
    print(f"Alpha error (with indecisive counting): {alpha_error:.8f}, Beta error (with indecisive counting): {beta_error:.8f}")
    alpha_error = sum(node['p1'] * node['count'] for node in leafs if node['decision'] == "p2")
    beta_error = sum(node['p2'] * node['count'] for node in leafs if node['decision'] == "p1")
    print(f"Alpha error (without indecisive counting): {alpha_error:.8f}, Beta error (without indecisive counting): {beta_error:.8f}")
# Example usage
# n = 5
# coinflips = 5  # Smaller for visualization purposes
# p1 = 0.5
# p2 = 0.6

# tree = create_tree(coinflips, p1, p2)

# fixed_tree = apply_fixed_sample_size_test(tree, coinflips = coinflips, bf_crit=1.33)
# print_table(fixed_tree)
# plot_tree_plotly(fixed_tree)
# opt_stop_tree = apply_optional_stopping(tree, bf_crit=1.33)
# print("\nAfter applying optional stopping:\n")
# print_tree(opt_stop_tree)
# print_table(opt_stop_tree)
# After applying optional stopping
# plot_tree_plotly(opt_stop_tree)

# Example
# n = 7
coinflips = 5
p1 = 0.5
p2 = 0.6
bf_crit = 1.33
#bf_crit2 = 1 / 1.05

tree = create_tree(coinflips, p1, p2)
#fig_tree = plot_tree_plotly(tree)
#fig_tree.write_image("figures/report/steele_replication_plain_tree.pdf")
other_tree = create_tree(coinflips, p1, p2)
fixed_tree = apply_fixed_sample_size_test(tree, coinflips = coinflips, bf_crit=bf_crit)
fixed_tree = apply_fixed_sample_size_test(tree, coinflips= 2, bf_crit=bf_crit)
print_table(tree)
fixed_fig = plot_tree_plotly(fixed_tree)
fixed_fig.write_image("figures/report/crossover_visualisation.pdf")
opt_stop_tree = apply_optional_stopping(other_tree, coinflips=coinflips, bf_crit=bf_crit)
print("\nAfter applying optional stopping:\n")
#print_tree(opt_stop_tree)
print_table(opt_stop_tree)
#fig = plot_tree_plotly(opt_stop_tree)

