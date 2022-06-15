import numpy as np


def dynamic_hamming(read: str, reference: str):
    """Scoring is callculated per each step the read walks on the reference

    Parameters
    ----------
    read : str
        meta-genomic read
    reference : str
        reference sequence

    Returns
    -------
    int, float
        score
    """

    start_idx = 0
    end_idx = len(read)
    scores = []
    if len(reference) > len(read):
        start_idx = 0
        end_idx = len(read)
        top_score = 0
        while end_idx - 1 < len(reference):
            ref_piece = reference[start_idx:end_idx]

            # score it
            score = hamming_distance_score(read=read, reference=ref_piece)
            if score == 0.0:
                scaled_score = 1.0 - score
                top_score = scaled_score
                break

            # increment
            start_idx += 1
            end_idx += 1
            scores.append(score)

    # if the read is larger than the actual gene, just swap
    else:
        print("Warning: A read is larger than the annotated gene, swapping...")
        # swapping
        small_gene = reference
        large_read = read

        # setting starting and ending window frame
        start_idx = 0
        end_idx = len(small_gene)

        top_score = 0
        while end_idx - 1 < len(large_read):
            ref_piece = large_read[start_idx:end_idx]

            # score it
            score = hamming_distance_score(read=small_gene, reference=ref_piece)
            if score == 0.0:
                print("Perfect match found")
                scores = score
                break

            # increment
            start_idx += 1
            end_idx += 1
            scores.append(score)

    if top_score == 1.0:
        return 1.0 - top_score
    elif isinstance(scores, float):
        return 1.0 - scores
    elif len(scores) > 0:
        return 1.0 - min(scores)


def hamming_distance_score(read: str, reference: str):
    """Calculates dissimilarity scores using hamming distance calculations

    Parameters
    ----------
    read : str
        meta-genomic read
    reference : str
        reference sequence

    Returns
    -------
    int, float
        hamming score
    """

    # converting into numpy arrays
    if isinstance(read, str):
        read = np.array([s for s in read])
    if isinstance(reference, str):
        reference = np.array([r for r in reference])

    read_ne_ref = read != reference
    score = np.average(read_ne_ref)
    return score
