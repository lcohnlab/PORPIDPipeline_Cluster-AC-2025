import pandas as pd
import yaml, sys, os, re

def create_config(input_csv_fhs):
    config = {}
    for input_csv in input_csv_fhs:
        dataset = {}
        #read in primer sheet
        df = pd.read_csv(input_csv)
        #filter to UMI samples
        df2 = df[(pd.notna(df["UMI"]))]
        for row in range(df2.shape[0]):
            if df2.iloc[row]["UMI"] == "dUMI":
                sec_str_primer = df2.iloc[row]["Second_Strand_Primer_Sequence"]
            else:
                sec_str_primer = df2.iloc[row]["Forward_Primer_2ndRd_Sequence"]
            
            if pd.isna(df2.iloc[row]["Panel"]):
                panel = "Error: No panel listed for this sample"
            else:
                panel = "panels/" + df2.iloc[row]["Panel"]    

            
            cDNA_primer = df2.iloc[row]["sUMI_Primer_Sequence"]
            m = re.search('N+', cDNA_primer)
            cDNA_primer = cDNA_primer[:m.start() - 6] + cDNA_primer[m.start()-6:m.start()].lower() +  cDNA_primer[m.start():]
									
            dataset[df2.iloc[row]["Sample"]] = {
                "cDNA_primer": cDNA_primer,
                "sec_str_primer": sec_str_primer,
                "panel": panel
                }

        config[os.path.basename(input_csv).split('.')[0]] = dataset
    sys.stdout.write(yaml.dump(config))

    filename = os.path.basename(input_csv).split('.')[0] + '-config.yaml'
    with open(filename, 'w') as file:
       doc = yaml.dump(config, file)
    sys.stdout.write('Moving yaml file to data/in/config directory\n')
    new_path = "../" + filename
    os.replace(filename, new_path)
    sys.stdout.write("Complete!\n")

create_config(sys.argv[1:])
