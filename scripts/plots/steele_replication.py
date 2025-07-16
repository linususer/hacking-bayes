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
        'decision': "indecisive" # Either "indecisive", "p1", or "p2"
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
              f"P1: {node['p1']:.4f}, P2: {node['p2']:.4f}, BF: {node['BF']}")
    print(f"Total nodes: {len(tree)}, Stopped nodes: {stopped_nodes}")

def apply_fixed_sample_size_test(tree, coinflips, bf_crit):
    tree_depth = coinflips
    for k in range(tree_depth + 1):
        idx = get_idx(tree_depth, k)
        node = tree[idx]
        bf = node['BF']
        if bf is not None:
            node['stopped'] = True  # Mark as "stopped" if BF crosses threshold
            if bf > bf_crit:
                node['decision'] = "p1"
            elif bf < 1 / bf_crit:
                node['decision'] = "p2"
            else:
                node['decision'] = "indecisive"
    return tree


def apply_optional_stopping(tree, bf_crit):
    for idx, node in tree.items():
        BF = node.get('BF')
        if BF is not None and (BF > bf_crit or BF < (1 / bf_crit)):
            parents = node['parent']
            p1_idx, p2_idx = parents
            if len(parents) == 1 or len(parents) == 2:
                # if parents are not stopped, we still count
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
                    node['count'] = tree[parents[0]]['count'] if parents[0] is not None else 0
                    if parents[1] is not None:
                        node['count'] = tree[parents[1]]['count']
            node['stopped'] = True
            if BF > bf_crit:
                node['decision'] = "p1"
            elif BF < (1 / bf_crit):
                node['decision'] = "p2"
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
                node['stopped'] = False
                node['decision'] = "indecisive"
    return tree

def plot_tree_plotly(tree):
    # Build edges and positions
    nodes = list(tree.keys())
    # labels for nodes with event values
    labels = [tree[n]['name'] for n in nodes]
    counts = [tree[n]['count'] for n in nodes]
    bfs = [tree[n]['BF'] for n in nodes]
    stopped = [tree[n]['stopped'] for n in nodes]
    
    # Calculate positions similarly as before for layout:
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
    
    # Build node scatter
    node_x = [pos[n][0] for n in nodes]
    node_y = [pos[n][1] for n in nodes]
    # Collect text annotations to the right of nodes
    text_x = []
    text_y = []
    text_labels = []
    node_color = []
    for n in nodes:
        x, y = pos[n]
        if tree[n]['count'] == 0:
            node_color.append('gray')
        elif tree[n]['stopped']:
            node_color.append('red')
            if tree[n]['decision'] == "p1":
                # Add independent text (right side of node)
                text_x.append(x + 0.15)      # Offset right
                text_y.append(y - 0.01)      # Slight Y adjustment
                text_labels.append(f"-")
            elif tree[n]['decision'] == "p2":
                # Add independent text (right side of node)
                text_x.append(x + 0.15)      # Offset right
                text_y.append(y - 0.01)      # Slight Y adjustment
                text_labels.append(f"+")
            else:
                text_x.append(x + 0.15)      # Offset right
                text_y.append(y - 0.01)      # Slight Y adjustment
                text_labels.append(f"?")
        else:
            node_color.append('skyblue')


    fig = go.Figure()

    # Add edges as lines
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y,
                             mode='lines',
                             line=dict(color='gray', width=1),
                             hoverinfo='none'
                            ))
    
    # Add nodes as markers
    fig.add_trace(go.Scatter(x=node_x, y=node_y,
                             mode='markers+text',
                             marker=dict(color=node_color, size=20, line_width=1),
                             text=labels,
                             textposition="top center",
                             hoverinfo='text',
                             hovertext=[f"Node {n}: {tree[n]['name']}<br>BF: {tree[n]['BF']}<br>Count: {tree[n]['count']}" for n in nodes]
                            ))
    
    # Add text annotations for counts
    fig.add_trace(go.Scatter(x=text_x, y=text_y,
                                mode='text',
                                text=text_labels,
                                textfont=dict(size=20, color='black'),
                                hoverinfo='none'
                                ))

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
    # add count number in the node
    for i, count in enumerate(counts):
        if count > 0:
            fig.add_annotation(
                x=node_x[i], y=node_y[i],
                text=str(count),
                showarrow=False,
                font=dict(size=14, color='black'),
                align='center',
            )
    fig.show()

# print table of tree nodes without 0 counts
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
                f"{node['p1']:.4f}\t{node['p2']:.4f}\t{round(node['BF'],3)}\t{node['decision']}")
    avg_tree_depth_p1 = sum([(node['event_value'][0] + node['event_value'][1]) * node['count'] * node['p1'] for node in tree.values()])
    avg_tree_depth_p2 = sum([(node['event_value'][0] + node['event_value'][1]) * node['count'] * node['p2'] for node in tree.values()])
    max_depth = max(node['event_value'][0] + node['event_value'][1] for node in tree.values())
    print(f"\nAverage tree depth for P(x | p1): {avg_tree_depth_p1:.4f} from max depth {max_depth}")
    print(f"\nAverage tree depth for P(x | p2): {avg_tree_depth_p2:.4f} from max depth {max_depth}")
    alpha_error = sum(node['p1'] * node['count'] for node in leafs if node['decision'] == "p2") + sum(node['p1'] * node['count'] for node in leafs if node['decision'] == "indecisive")
    beta_error = sum(node['p2'] * node['count'] for node in leafs if node['decision'] == "p1") + sum(node['p2'] * node['count'] for node in leafs if node['decision'] == "indecisive")
    print(f"Alpha error (with indecisive counting): {alpha_error:.4f}, Beta error (with indecisive counting): {beta_error:.4f}")
    alpha_error = sum(node['p1'] * node['count'] for node in leafs if node['decision'] == "p2")
    beta_error = sum(node['p2'] * node['count'] for node in leafs if node['decision'] == "p1")
    print(f"Alpha error (without indecisive counting): {alpha_error:.4f}, Beta error (without indecisive counting): {beta_error:.4f}")
# Example usage
# n = 5
coinflips = 5  # Smaller for visualization purposes
p1 = 0.5
p2 = 0.6

tree = create_tree(coinflips, p1, p2)

# Usage with the previously generated tree:
fixed_tree = apply_fixed_sample_size_test(tree, coinflips = coinflips, bf_crit=1.33)
print_table(fixed_tree)
# plot_tree_plotly(fixed_tree)
# opt_stop_tree = apply_optional_stopping(tree, bf_crit=1.33)
# print("\nAfter applying optional stopping:\n")
# print_tree(opt_stop_tree)
# print_table(opt_stop_tree)
# After applying optional stopping
# plot_tree_plotly(opt_stop_tree)

# Example
# n = 7
coinflips = 7
p1 = 0.5
p2 = 0.6
bf_crit = 1.33

tree = create_tree(coinflips, p1, p2)
fixed_tree = apply_fixed_sample_size_test(tree, coinflips = coinflips, bf_crit=bf_crit)
#print_table(tree)
#plot_tree_plotly(fixed_tree)
opt_stop_tree = apply_optional_stopping(tree, bf_crit=bf_crit)
print("\nAfter applying optional stopping:\n")
#print_tree(opt_stop_tree)
print_table(opt_stop_tree)
#plot_tree_plotly(opt_stop_tree)

