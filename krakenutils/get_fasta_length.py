import subprocess

def get_fasta_length(fasta: str) -> int | None:
    """
    Returns the number of fasta entries in a file by counting the character ">"

    Args:
        fasta (string): Path to a fasta file

    Returns:
        integer: Count of ">" characters in fasta
    """
    cmd = (["grep", "-c", ">" , fasta])
    sp = subprocess.run(cmd, capture_output=True, text=True)
    if sp.returncode == 0:
        return int(sp.stdout.strip("\n"))
    else:
        print(sp.stderr)
        exit()