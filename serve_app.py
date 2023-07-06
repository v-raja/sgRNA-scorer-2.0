import gradio as gr
from identifyAndScore import runPipeline, validPAM
from ray import serve
from ray.serve.gradio_integrations import GradioIngress, GradioServer

# import pydevd_pycharm
# pydevd_pycharm.settrace('localhost', port=121, stdoutToServer=True, stderrToServer=True)
# Create a builder function for the Gradio Interface
def build_gradio_app():
    def run_scorer(input_file, output_file, pam_orientation, spacer_length, pam_sequence):
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
        runPipeline(input_file, output_file, spacer_length, pam_sequence, pam_orientation)

        # Return output
        with open(output_file, 'r') as file:
            return file.read()

    iface = gr.Interface(
        fn=run_scorer,
        inputs=[
            gr.inputs.File(label="Input File"),
            gr.inputs.Textbox(label="Output File"),
            gr.inputs.Dropdown(choices=["5", "3"], label="PAM Orientation"),
            gr.inputs.Textbox(label="Spacer length"),
            gr.inputs.Textbox(label="PAM Sequence")
        ],
        outputs="text",
        examples=[
            ["eGFP.fasta", "eGFP.SpCas9.tab", "3", "20", "NGG"],
            ["eGFP.fasta", "eGFP.SaCas9.tab", "3", "21", "NNGRRT"],
            ["eGFP.fasta", "eGFP.Cpf1.tab", "5", "20", "TTTN"]
        ]
    )

    return iface

    # iface.add_selection_field("Preset", options=presets.keys(), type="run_scorer", action=select_preset)

# Generic Gradio server
entry_point_for_this_gradio_app = GradioServer.options(ray_actor_options={"num_cpus": 1, "num_gpus": 0}).bind(
    build_gradio_app
)