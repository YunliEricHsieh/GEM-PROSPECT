import subprocess
import os

# --- CONFIG ---
input_csv = 'Results/CataPro/enzymes_and_substrates.csv'
output_kcat = 'Results/CataPro/enzymes_and_substrates_kcats.csv'

# Verification: Ensure the input file exists
if not os.path.exists(input_csv):
    print(f"🔴 Error: {input_csv} not found. Did you run the preparation script?")
else:
    print(f"🚀 Starting CataPro inference on curated data...")
    
    catapro_cmd = [
        "python", "CataPro/inference/predict.py",
        "-inp_fpath", input_csv,
        "-model_dpath", "CataPro/models",
        "-batch_size", "128",
        "-device", "cpu",
        "-out_fpath", output_kcat
    ]

    try:
        subprocess.run(catapro_cmd, check=True)
        print(f"\n✅ Done! Predicted kcats saved to: {output_kcat}")
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Prediction failed: {e}")