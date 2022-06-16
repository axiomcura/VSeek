from pathlib import Path
from vseek.common.errors import FastaFileNotFound
from vseek.utils.data_structs import ReadRecord


class SequenceIO:
    def __init__(self, sequence_path):
        # checking
        check = Path(sequence_path).is_file()
        if check is False:
            raise FastaFileNotFound(f"Unable to find {sequence_path}")

        self.sequence_path = sequence_path

    def lazy_load_fasta(self):
        """Loads all fasta sequences as FastaReadRecords

        Yields
        ------
        _type_
            _description_
        """
        with open(self.sequence_path, "r") as fasta_file:
            entry = []
            for idx, line in enumerate(fasta_file):

                cleaned_line = line.strip().replace("\n", "")

                # append the first line only
                if idx == 0 and cleaned_line.startswith(">"):
                    entry.append(cleaned_line)
                    continue

                elif idx != 0 and cleaned_line.startswith(">"):
                    # convert the chunk into a single
                    fasta_record = self._convert(entry)

                    # reassign list entry as empty
                    entry = []

                    yield fasta_record

                entry.append(cleaned_line)

    def _convert(self, entry: list[str]):
        """converts entry into a FastaReadRecord

        Parameters
        ----------
        entry : list[str]
            list containing header and sequence

        Returns
        -------
        ReadRecord
            Datatype that contains header_id, fragment_id, sequence and its length
        """
        entry_data = tuple(entry[0].split())
        header_id = entry_data[0]
        frag_id = entry_data[1]
        sequence = "".join(entry[1:])
        length = len(sequence)

        return ReadRecord(
            srr_id=header_id, fragment_id=frag_id, sequence=sequence, length=length
        )
