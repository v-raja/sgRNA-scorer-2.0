import gradio as gr
from identifyAndScore import runPipeline, validPAM  # Assuming this code is saved in sgRNA_Scorer.py


def run_scorer(input_path, output_path, spacer_length, pam_sequence, pam_orientation):
    # Check PAM validity
    if not validPAM(pam_sequence):
        return 'PAM sequence has invalid characters. PAM sequence must only contain A,C,T,G,K,M,R,Y,S,W,B,V,H,D,N'

    # Check spacer length
    if not spacer_length.isdigit() or int(spacer_length) < 14:
        return 'Spacer length must be an integer with a minimum size of 14'

    # Check PAM orientation
    if pam_orientation not in ['5', '3']:
        return 'Valid PAM orientations are 5 or 3'

    # Perform operation
    runPipeline(input_path, output_path, spacer_length, pam_sequence, pam_orientation)

    # Return output
    with open(output_path, 'r') as output_file:
        return output_file.read()


iface = gr.Interface(
    fn=run_scorer,
    inputs=["text", "text", "text", "text", "text"],
    outputs="text")

iface.launch()
