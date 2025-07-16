# create a pascal-triangle-like tree based on the number of coinflips
# for probability p1 and p2
# indexing with idx(tree_depth,k) = (tree_depth+1)tree_depth/2 + k + 1
# where k is the "column"

createTree <- function(coinflips, p1, p2) {
    tree <- list()
    tree_depth <- 0
    idx <- getIdx(tree_depth, 0)
    tree <- addNode(tree, idx, c(NULL, NULL), c(0, 0), p1, p2)
    tree_depth <- 1
    # iterate for every coinflip
    for (i in 1:coinflips) {
            # iterate for every layer (nodes = tree_depth + 1)
            for (k in 0:tree_depth) {
                idx <- getIdx(tree_depth, k)
                parents <- getParents(tree_depth, k)
                event_seq <- c(tree_depth - k, k)
                tree <- addNode(tree, idx, parents, event_seq, p1, p2)
            }
            # increment tree depth
            tree_depth <- tree_depth + 1
        }
    return(tree)
    }


# get the parent nodes given the indices of the node
getParents <- function(tree_depth, k) {
    # if node has a right neighbour
    if ((k + 1) <= tree_depth) {
        parent1 <- getIdx(tree_depth - 1, k)
    } else { # if node is most right
        parent1 <- NULL
    }
    # if node has a left neighbour
    if ((k - 1) >= 0) {
        parent2 <- getIdx(tree_depth - 1, k - 1)
    } else { # if node is most left
        parent2 <- NULL
    }
    c(parent1, parent2)
}
# get the index of the node
getIdx <- function(tree_depth, k) {
    (((tree_depth + 1) * tree_depth) / 2) + k + 1
}
# add node to tree at idx, parent nodes, event_seq (HEADS,TAILS), p1, p2
addNode <- function(tree, idx, parents, event_seq, p1, p2) {
        count <- choose(event_seq[1] + event_seq[2], event_seq[1]) # binomial coefficient
        p1 <- count * p1^event_seq[1] * (1 - p1)^event_seq[2] # probability for p1
        p2 <- count * p2^event_seq[1] * (1 - p2)^event_seq[2]  # probability for p2
        BF <- p1 / p2 # Bayes Factor
        node <- list(
        name = paste0("H", event_seq[1], "T", event_seq[2]),
        parent = parents, # c(parents[1], parents[2]),
        event_value = event_seq, # c(HEADS, TAILS),
        count = count,
        p1 = p1,
        p2 = p2,
        BF = BF,
        stopped = FALSE
    )
    tree[[idx]] <- node
    return(tree)
}


printTree <- function(tree) {
    stopped_nodes <- 0
    for (i in seq_along(tree)) {
        node <- tree[[i]]
        if (node$stopped) {
            stopped_nodes <- stopped_nodes + 1
        }
        cat(sprintf("Node %d: %s, Parents: %s, Event Value: %s, Count: %.0f, P1: %.4f, P2: %.4f, BF: %.4f\n",
                    i, node$name, paste(node$parent, collapse = ", "),
                    paste(node$event_value, collapse = ", "), node$count,
                    node$p1, node$p2, node$BF))
    }
    # print number of nodes with stopped = TRUE
    cat(sprintf("Total nodes: %d, Stopped nodes: %d\n", length(tree), stopped_nodes))
}
tree <- createTree(coinflips = 50, p1 = 0.5, p2 = 0.6)
printTree(tree)

# Given a tree and a critical bayes factor threshold
# apply the Optional Stopping Rule onto the tree
applyOptionalStopping <- function(tree, bf_crit){
    for (idx in seq_along(tree)) {
        node <- tree[[idx]]
        
        # Check if node BF is outside threshold
        if (!is.null(node$BF) && (node$BF > bf_crit || node$BF < 1 / bf_crit)) {
            node$stopped <- TRUE
            node$count <- 0
            node$BF <- NULL
        } else {
            parents <- node$parent
            if (!is.null(parents) && length(parents) == 2) {
                p1 <- parents[1]
                p2 <- parents[2]
                
                p1_stopped <- !is.null(tree[[p1]]) && tree[[p1]]$stopped
                p2_stopped <- !is.null(tree[[p2]]) && tree[[p2]]$stopped
                
                if (p1_stopped && p2_stopped) {
                    node$stopped <- TRUE
                    node$count <- 0
                    node$BF <- NULL
                } else {
                    node$p1 <- node$p1 / node$count
                    node$p2 <- node$p2 / node$count
                    node$count <- tree[[p1]]$count + tree[[p2]]$count
                    node$p1 <- node$count * node$p1
                    node$p2 <- node$count * node$p2
                    node$BF <- node$p1 / node$p2
                }
            }
        }
        tree[[idx]] <- node
    }

    return(tree)
}


opt_stop_tree <- applyOptionalStopping(tree, bf_crit = 3)
printTree(opt_stop_tree)
