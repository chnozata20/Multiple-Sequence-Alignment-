class Color:
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    END = '\033[0m'
class TreeNode:
    def __init__(self, data):
        self.data = data
        self.children = []
        self.seq = []

    def add_child(self, child):
        self.children.append(child)

def depth_first_search(node, target):
    if node.data == target:
        return node

    for child in node.children:
        result = depth_first_search(child, target)
        if result:
            return result
def print_tree(node, prefix="", is_last=True):
    print(prefix + ("└── " if is_last else "├── ") + node.data)
    child_count = len(node.children)
    for i, child in enumerate(node.children):
        is_last_child = i == child_count - 1
        child_prefix = prefix + ("    " if is_last else "│   ")
        print_tree(child, child_prefix, is_last_child)
def breadth_first_search(root, target):
    queue = [root]
    while queue:
        node = queue.pop(0)
        if node.data == target:
            return node
        queue.extend(node.children)

def postfix_traversal(node):
    result = []
    if node is not None:
        for child in reversed(node.children):
            result.extend(postfix_traversal(child))
        result.append(node.data)

    return result

def add_seq_to_tree(node,sequences,counter,blosum62_matrix,gap_penalty):
    if node is not None:
        for child in reversed(node.children):
            if len(child.data.split("*")) < 2:
                for index in range(0, len(sequences), 2):
                    if sequences[index] == child.data:
                        child.seq=[sequences[index+1]]
            add_seq_to_tree(child,sequences,counter,blosum62_matrix,gap_penalty)
        if len(node.data.split("*")) == counter:
            node.seq=(multiple_alignment(node.children[0].seq,node.children[1].seq, blosum62_matrix, gap_penalty=gap_penalty ))




def global_alignment(seq1, seq2, blosum62_matrix, gap_penalty=-5):
    # Initialize the alignment matrix
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    alignment_matrix = [[0] * cols for _ in range(rows)]

    # Fill the alignment matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match = alignment_matrix[i-1][j-1] + blosum62_matrix[seq1[i-1]][seq2[j-1]]
            delete = alignment_matrix[i-1][j] + gap_penalty
            insert = alignment_matrix[i][j-1] + gap_penalty
            alignment_matrix[i][j] = max(match, delete, insert)

    # Traceback to find the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = rows - 1, cols - 1

    while i > 0 or j > 0:
        if i > 0 and alignment_matrix[i][j] == alignment_matrix[i-1][j] + gap_penalty:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "*" + aligned_seq2
            i -= 1
        elif j > 0 and alignment_matrix[i][j] == alignment_matrix[i][j-1] + gap_penalty:
            aligned_seq1 = "*" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
    return aligned_seq1, aligned_seq2

def multiple_alignment(seq2, seq1, blosum62_matrix, gap_penalty=-5):
    # Initialize the alignment matrix

    rows = len(seq1[0]) + 1
    cols = len(seq2[0]) + 1
    alignment_matrix = [[0] * cols for _ in range(rows)]
    alignment_matrix_for_goal = [[0] * cols for _ in range(rows)]

    # Fill the alignment matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match_value = 0
            delete_value = 0
            insert_value = 0
            for k in range(0,len(seq1)):
                for n in range(0, len(seq2)):
                    if seq1[k][i-1] != "*" and seq2[n][j-1] != "*":
                        match_value += blosum62_matrix[seq1[k][i-1]][seq2[n][j-1]]
                        delete_value += gap_penalty
                        insert_value += gap_penalty
                        if len(seq1) == 1 and len(seq2) == 1:
                            delete_value = gap_penalty
                            insert_value = gap_penalty
                    elif seq1[k][i-1] == "*" and seq2[n][j-1] == "*":
                        match_value += blosum62_matrix[seq1[k][i-1]][seq2[n][j-1]]
                        delete_value += blosum62_matrix[seq1[k][i-1]][seq2[n][j-1]]
                        insert_value += blosum62_matrix[seq1[k][i-1]][seq2[n][j-1]]
                    else:
                        match_value += gap_penalty
                        if seq1[k][i-1] != "*":
                            delete_value += gap_penalty
                        else:
                            delete_value = blosum62_matrix[seq1[k][i - 1]]["*"]
                        if seq2[n][j-1] != "*":
                            insert_value += gap_penalty
                        else:
                            insert_value = blosum62_matrix["*"][seq2[n][j-1]]

            max_temp = max(alignment_matrix[i-1][j-1] + match_value,
                           alignment_matrix[i-1][j] + delete_value,
                           alignment_matrix[i][j-1] + insert_value)

            alignment_matrix[i][j] = max_temp
            if max_temp == alignment_matrix[i-1][j-1] + match_value:
                max_temp_2 = match_value
            elif max_temp == alignment_matrix[i-1][j] + delete_value:
                max_temp_2 = delete_value
            else:
                max_temp_2 = insert_value

            alignment_matrix_for_goal[i][j] = max_temp_2

    # Traceback to find the aligned sequences
    aligned_seq1 = []
    aligned_seq2 = []
    for x in range(len(seq1)):
        aligned_seq1.append("")
    for x in range(len(seq2)):
        aligned_seq2.append("")

    i, j = rows - 1, cols - 1

    while i > 0 or j > 0:
        if i > 0 and j > 0 and alignment_matrix[i][j] == alignment_matrix[i-1][j-1] + alignment_matrix_for_goal[i][j]:
            for index in range(len(seq1)):
                aligned_seq1[index] = seq1[index][i - 1] + aligned_seq1[index]
            for index in range(len(seq2)):
                aligned_seq2[index] = seq2[index][j - 1] + aligned_seq2[index]
            i -= 1
            j -= 1
        elif i > 0 and alignment_matrix[i][j] == alignment_matrix[i-1][j] + alignment_matrix_for_goal[i][j]:
            for index in range(len(seq1)):
                aligned_seq1[index] = seq1[index][i-1] + aligned_seq1[index]
            for index in range(len(seq2)):
                aligned_seq2[index] = "*" + aligned_seq2[index]
            i -= 1
        elif j > 0 and alignment_matrix[i][j] == alignment_matrix[i][j-1] + alignment_matrix_for_goal[i][j]:
            for index in range(len(seq1)):
                aligned_seq1[index] = "*" + aligned_seq1[index]
            for index in range(len(seq2)):
                aligned_seq2[index] = seq2[index][j - 1] + aligned_seq2[index]
            j -= 1
    return aligned_seq1+aligned_seq2

def read_blosum62_matrix(file_path="./txt/Blosum62.txt"):
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Parse the BLOSUM62 matrix
    amino_acids = lines[0].split()
    blosum62_matrix = {}
    for line in lines[1:]:
        values = line.split()
        row_amino_acid = values[0]
        row_scores = [int(score) for score in values[1:]]
        blosum62_matrix[row_amino_acid] = dict(zip(amino_acids, row_scores))

    return blosum62_matrix


def compute_similarity_score(aligned_seq1, aligned_seq2):
    # Calculate the number of exact matches
    exact_matches = sum(a == b and a != "*" for a, b in zip(aligned_seq1, aligned_seq2))

    print("exact matches: ", exact_matches)
    # Calculate the aligned sequence length
    if len(aligned_seq1) >= len(aligned_seq2):
        aligned_length = len(aligned_seq1)
    else:
        aligned_length = len(aligned_seq2)
    print("aligned length: ", aligned_length)
    # Calculate the similarity score
    similarity_score = exact_matches / aligned_length

    return similarity_score


def main():
    # Get gap penalty from the user
    gap_penalty = int(input("Enter the gap penalty: "))

    # Read BLOSUM62 matrix
    blosum62_matrix = read_blosum62_matrix()
    # Read sequences from an input file
    input_file = "./txt/Input.txt"  # Provide the path to your input file
    with open(input_file, "r") as file:
        lines = file.readlines()
        sequences = [line.strip().upper().replace(">", "") for line in lines]

    k = int(len(sequences)/2)
    similarity_matrix = [[0.0] * k for _ in range(k)]
    asd = {}
    aligned_seq1 = ""
    aligned_seq2 = ""

    # Perform pairwise sequence alignments
    for i in range(1, len(sequences), 2):
        for j in range(i + 2, len(sequences), 2):
            sequence1 = sequences[i]
            sequence2 = sequences[j]
            aligned_seq1, aligned_seq2 = global_alignment(sequence1, sequence2, blosum62_matrix, gap_penalty=gap_penalty)
            print(f"{Color.RED}Before Pairwise Alignment:{Color.END}")
            print(int(i/2)+1, "(v) and (w)", int(j/2)+1, " : ","\n", sequence1,"\n", sequence2, sep="")
            print(f"{Color.YELLOW}After Pairwise Alignment:{Color.END}")
            print("v:", aligned_seq1)
            print("w:", aligned_seq2)
            # Display the aligned sequences for each pair
            print(f"Alignment for sequences {int(i/2)+1} and {int(j/2)+1}:")
            similarity_score = compute_similarity_score(aligned_seq1, aligned_seq2)
            similarity_matrix[int(i/2)][int(j/2)] = similarity_score
            similarity_matrix[int(j/2)][int(i/2)] = similarity_score  # Since it's a symmetric matrix
            print("\n" + "=" * 30 + "\n")  # Separate alignments for clarity

    # Perform similarity calculation
    for i in range(0, len(sequences), 2):
        for j in range(i + 2, len(sequences), 2):
            seq1 = sequences[i + 1]
            seq2 = sequences[j + 1]
            print(int(i/2)+1, " ve ", int(j/2)+1, " : ", seq1, seq2)
            # similarity_score = compute_similarity_score(aligned_seq1, aligned_seq2)
            # similarity_matrix[int(i/2)][int(j/2)] = similarity_score
            # similarity_matrix[int(j/2)][int(i/2)] = similarity_score  # Since it's a symmetric matrix

            key = (sequences[i], sequences[j])
            # asd[key] = similarity_score

    # Display the similarity matrix
    # Display the similarity matrix
    print()
    print("Similarity Matrix:")
    max_seq_len = max(len(seq) for seq in sequences[::2])  # En uzun başlık uzunluğunu bul
    header = " " * max_seq_len + " " .join(seq.rjust(max_seq_len) for seq in sequences[::2])  # Başlık için
    print(header)
    for i, row in enumerate(similarity_matrix):
        formatted_row = [f"{value:.3f}" for value in row]
        print(f"{sequences[i * 2].rjust(max_seq_len)} {' '.join(formatted_row)}")

    """
    print("Similarity Matrix:")
    for key, value in asd.items():
        print(f"{key}: {value}")
    """

    labels = []
    for i in range(0, len(sequences), 2):
        labels.append(sequences[i])




    def UPGMA(similarity_matrix, labels):

        def find_min(matrix):
            path = ""
            max_key = ""
            max_key_2 = ""
            max = float('-inf')

            for key, value in matrix.items():
                for key_2, value_2 in matrix[key].items():
                    if key != key_2 and value_2 > max:
                        max = value_2
                        path=str(key_2)+"*"+str(key)
                        max_key=key
                        max_key_2=key_2


            temp = {}
            for key, value in matrix[max_key].items():
                temp[key] = (float(matrix[max_key][key])+float(matrix[max_key_2][key]))/2

            matrix[path] = temp

            for key, value in matrix.items():
                if key == path:
                    matrix[key][path] = 0
                else:
                    matrix[key][path] = matrix[path][key]

            keys_to_delete = []

            for key, value in matrix.items():
                if key == max_key or key == max_key_2:
                    keys_to_delete.append(key)
                else:
                    keys_to_delete_2 = []
                    for key_2, value_2 in matrix[key].items():
                        if key_2 == max_key or key_2 == max_key_2:
                            keys_to_delete_2.append(key_2)

                    for key_2 in keys_to_delete_2:
                        del matrix[key][key_2]

            for key in keys_to_delete:
                del matrix[key]

            if len(matrix)>1:
                tree = find_min(matrix)

                node = breadth_first_search(tree, path)

                node.add_child(TreeNode(max_key))
                node.add_child(TreeNode(max_key_2))
                return tree
            else:
                tree = TreeNode(path)
                tree.add_child(TreeNode(max_key))
                tree.add_child(TreeNode(max_key_2))
                return tree


        #print(similarity_matrix)
        dict_matrix = {}
        for index in range(len(labels)):
            temp = {}
            for index_2 in range(len(similarity_matrix[index])):
                temp[labels[index_2]] = similarity_matrix[index][index_2]

            dict_matrix[labels[index]] = temp

        tree = find_min(dict_matrix)
        print()
        print("UPGMA Tree")
        print(print_tree(tree))

        for counter in range(2,(len(tree.data.split("*"))+1)):
            add_seq_to_tree(tree,sequences,counter,blosum62_matrix,gap_penalty)

        print()
        print(f"{Color.RED}Multiple Alignment Result{Color.END}")
        label_list = tree.data.split("*")
        for index in range(len(tree.seq)):
            print(f"{Color.YELLOW}{label_list[index]}{Color.YELLOW}")
            print(f"{Color.BLUE}{tree.seq[index]}{Color.BLUE}")

        postfix_result = postfix_traversal(tree)
        #print(postfix_result)
    UPGMA(similarity_matrix, labels)


main()