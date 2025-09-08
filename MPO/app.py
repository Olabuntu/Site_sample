import streamlit as st
from streamlit_option_menu import option_menu
import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
import os
import numpy as np
import pickle
import sklearn


def delete_all_files():
    for file in os.listdir():
        if file.endswith('.xlsx') or file.endswith('.csv'):
            os.remove(file)



def check(a):
    signs = ['<', '>', '<=', '>=']
    for item in a:
        for key, value in item.items():
        
            if value['lower'] is not None and not isinstance(value['lower'], (int, float)):
                return False
    
            if (not isinstance(value['upper'], (int, float)) and 
                value['upper'] not in signs):
                return False
            
                
            else:
                return True

def lipi(i):
    MW, HBA, HBD,  LogP, scores,Lipinski = 0,0,0,0,0,None
    score = 0
  
    mol = Chem.MolFromSmiles(i)
    if mol is None:
        MW, HBA, HBD, LogP, scores, Lipinski = None, None, None, None, None, "Invalid SMILES"
    else:
        MW = Descriptors.MolWt(mol)
        if Descriptors.MolWt(mol) > 500:
            score += 1
        HBA= Descriptors.NumHAcceptors(mol)
        if HBA > 10:
            score += 1
        HBD = Descriptors.NumHDonors(mol)
        if HBD > 5:
            score += 1
        LogP = Descriptors.MolLogP(mol)
        if LogP > 5:
            score += 1
        scores = score
        if score > 1:
            Lipinski = 'fail'
        else:
            Lipinski = 'pass'


    return MW, HBA, HBD, LogP, scores, Lipinski

def Lipinski(smile):
    MW, HBA, HBD, LogP, scores,Lipinski = [], [], [], [], [], []
    if type(smile) is str:
        MW, HBA, HBD, LogP, scores, Lipinski = lipi(smile)
        return MW, HBA, HBD, LogP, scores, Lipinski


    elif type(smile) is list:
        for i in smile:
            MW_i, HBA_i, HBD_i, LogP_i, scores_i, Lipinski_i = lipi(i)
            MW.append(MW_i)
            HBA.append(HBA_i)
            HBD.append(HBD_i)
            LogP.append(LogP_i)
            scores.append(scores_i)
            Lipinski.append(Lipinski_i)


        ids = smile
        df = pd.DataFrame({'smiles':ids, 
                           'Molecular Weight': MW,
                           'HBA': HBA,
                           'HBD': HBD,
                           'LogP': LogP,
                           'Failed Criteria': scores,
                           'Comment': Lipinski})
        return df


    
    
def file_formation (file_data):
    # Create a DataFrame from the input data
   
    file_data = [col for col in file_data]
    file_data.insert(0, 'ID')
    df = pd.DataFrame(columns= file_data)

    # Save the DataFrame to a CSV file
    csv_file = "input_data.csv"
    df.to_csv(csv_file, index=False)

    return True


def analysis():
        label = []
        limits = []
        weights = []
        file_name = 'uploaded_file.csv'
        try:
            df = pd.read_csv(file_name)
            ids =  df['ID'].tolist() 
            df_input =  df.copy()
            df.pop('ID')
            
            column = [col.lower()  for col in df.columns.tolist()]
            
            df.columns = column
            input_data1 = st.session_state.get("input_data", {})
            input_data = {col.lower(): input_data1[col] for col in input_data1.keys() if col.lower() in df.columns}
            keys = [col.lower() for col in  list(input_data.keys())]
        
        

            
            for key in keys:

                if 'range' in input_data[key]:
                    
                    lower = input_data[key]['range']['lower']
                    upper = input_data[key]['range']['upper']
                    weight = input_data[key]['range']['weight']
                    if upper == '<':
                        df[key] = df[key].where(df[key] < lower, 0)
                    elif upper == '>':
                        df[key] = df[key].where(df[key] > lower, 0)
                    elif upper == '<=':
                        df[key] = df[key].where(df[key] <= lower, 0)
                    elif upper == '>=':
                        df[key] = df[key].where(df[key] >= lower, 0)
                    #min max normalization
                    df[key] = df[key].fillna(0)
                    if df[key].max() != df[key].min():
                        df[key] = (df[key] - df[key].min()) / (df[key].max() - df[key].min())

                    df[key] = df[key] * weight
                        
                elif 'limits' in input_data[key]:
                    lower = input_data[key]['limits']['lower']
                    upper = input_data[key]['limits']['upper']
                    weight = input_data[key]['limits']['weight']
                    df[key] = df[key].where((df[key] >= lower) & (df[key] <= upper), 0)

                    df[key] = df[key].fillna(0)
                    if df[key].max() != df[key].min():
                        df[key] = (df[key] - df[key].min()) / (df[key].max() - df[key].min())

                    df[key] = df[key] * weight
            df.columns = list(input_data1.keys())
            df['MPO_score'] = df.sum(axis=1)
            df.insert(0,'ID',ids)
            df.sort_values(by='MPO_score', ascending =False, inplace=True)
            df.reset_index(drop=True, inplace=True)

            
            
            for y,x in  zip(input_data1.keys(),input_data1.values()):
        
                if 'range' in x.keys():
                    num = x['range']['lower']
                    sign = x['range']['upper']
                    limit = f'{sign} {num}'
                    weight = x['range']['weight']
                    
                
                elif 'limits' in x.keys():
                    lower = x['limits']['lower']
                    upper = x['limits']['upper']
                    limit = f'lower = {lower} and upper = {upper}'
                    weight = x['limits']['weight']

                label.append(y)
                limits.append(limit)
                weights.append(weight) 

            df2 = pd.DataFrame({'Parameters':label, 'Limits':limits,'Weights': weights})
            

        
            with pd.ExcelWriter('analysis_results.xlsx') as writer:
                df_input.to_excel(writer, sheet_name='Input_Data', index=False)
                df.to_excel(writer, sheet_name='MPO_Score', index=False)
                df2.to_excel(writer, sheet_name='Limits', index=False)


        except Exception as e:
            st.error(f"An error occurred during analysis, please check your input file. Error: {e}")

        



def check_columns(file_data, expected_columns):
    
    file_data = [col.lower() for col in file_data.columns]
    if file_data == expected_columns:
        
        return True 
    else: 
        return False
     

def descriptors(df, model):
    
    embeddings = []
    for smile in df['canonical_smiles']:
        molecule = Chem.MolFromSmiles(smile)
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, radius=2, nBits=2048)
        fingerprint_array = np.array(fingerprint)    
        embeddings.append(fingerprint_array)
    df1 = pd.DataFrame(embeddings)
    model = pickle.load(open(model, 'rb'))
    predictions = model.predict(df1)
    df2 = df.id
    df3 = pd.concat([df2, pd.Series(predictions, name='predictions')], axis=1)
    df3.columns = ['ID', 'Predicted pIC₅₀']
    df3['Predicted IC₅₀(nM)'] = 10 ** (-df3['Predicted pIC₅₀'])
    df3['Predicted IC₅₀']  = df3['Predicted IC₅₀']  * (10**9)
    return df3







st.set_page_config(
    page_title="Streamlit App",
    page_icon="🔍",
     layout="centered",
    initial_sidebar_state="collapsed"
)


# Sidebar menu
with st.sidebar:
    selected = option_menu(
        menu_title="Main Menu",
        options=["MPO", "Lipinski","Models","About", "Contact"],
        icons=["house", "info-circle", "envelope"],
        default_index=0,
        orientation="vertical",
        styles={
            "container": {"padding": "0!important"},
            "icon": {"color": "#ffffff"},
            "nav-link": {
                "font-size": "1.2rem",
                "text-align": "left",
                "margin": "0px",
                "--hover-color": "#f0f0f0"
            },
            "nav-link-selected": {"background-color": "#007bff"}
        }
    )

if "last_selected" not in st.session_state:
    st.session_state.last_selected = selected

if st.session_state.last_selected != selected:
    # Reset whatever needs clearing
    for key in list(st.session_state.keys()):
        if key not in ["last_selected"]:
            del st.session_state[key]
    st.session_state.last_selected = selected



if selected == "MPO":
    
    if "home" not in st.session_state:
        st.session_state["home"] = 1

    # Step 2: Show Parameter Selection
    if st.session_state.get("home") == 1:
        delete_all_files()
        parameters = [
            "Experimental IC₅₀", "Experimental KI", "Predicted pIC₅₀", "Predicted IC₅₀", "Predicted KI", "Predicted pKI", "FEP ΔΔG", "MM‑GBSA ΔG",
            "Docking Score", "Lipophilic Efficiency (LiPE)", "LogP", "LogS (Solubility)",
            "Molecular Weight", "TPSA", "Rotatable Bonds", "BBB Permeability",
            "Plasma Protein Binding", "CYP450 Inhibition", "hERG Inhibition",
            "AMES Toxicity", "PAINS & Structural Alerts", "Hydrogen Bond Donors", "Hydrogen Bond Acceptors",  
            "Synthetic Accessibility Score"
        ]

        st.title("Select all Parameters")
        options = st.multiselect(
            "Click to select Parameters to Display",
            options=parameters,
            key="parameter_selection"
        )

        if st.button("Submit", key="submit_button", help="Click to submit your parameter selections"):
            if not options:
                st.warning("⚠️ Please select at least one parameter to display.")
            else:
                st.session_state["selected_parameters"] = options
                st.session_state["submitted"] = True
                st.session_state["confirmed"] = False

        if st.session_state.get("submitted"):
            st.session_state["home"] = 2
            st.rerun()

    # Step 3: After submission, show confirmation options
    elif st.session_state.get("home") == 2:
        options = st.session_state.get("selected_parameters", [])
        st.title("Selected Parameters")
        st.success(
            f"You have selected {len(options)} parameters:\n\n" +
            "\n".join(f"- {opt}" for opt in options)
        )

        col1, col2, col3 = st.columns(3)

        with col1:
            st.button('🔄 Reset', key='reset_button', on_click=lambda: st.session_state.clear())

        with col3:
            if st.button('✅ Confirm', key='confirm_button'):
                with st.spinner("Loading..."):
                    st.session_state["confirmed"] = True
                    st.session_state["home"] = 3
                    st.rerun()

    # Step 4: After confirmation, go to Next phase
    elif st.session_state.get("confirmed") and st.session_state.get("home") == 3:
        input_data = {}
        weights = []
        total_weight = 0.0
        comparison_options = {
            "Greater than": ">",
            "Less than": "<",
            "Greater than or equal to": ">=",
            "Less than or equal to": "<="
        }
            
        st.markdown(
            """
            <style>
            .centered-title {
            text-align: center;
            font-size: 2rem;
            font-weight: bold;
            margin-top: 20px;
            margin-bottom: 20px;
            }
            </style>
            <div class="centered-title">Parameter Input Example</div>
            """,
            unsafe_allow_html=True
        )
        col1, col2,col3 = st.columns(3)
        with col1:
            st.markdown("<strong><b>Range</b></strong> eg: Greater than 1000", unsafe_allow_html=True)
        with col3:
            st.markdown("Limits eg: Lower limit 1000, Upper limit 2000")

        
        for parameter in st.session_state.get("selected_parameters", []):
            
            st.markdown(f"**{parameter}**")

            types = st.radio(
                "Select parameter type",
                options=["Range", "Limits"],
                key=f"type_{parameter}",
                horizontal=True)
            col1, col2, col3 = st.columns(3)
            
            if types == "Range":
                with col1:
                    dropdown = st.selectbox(
                        f'Select {parameter} Range',
                        options=["Greater than", "Less than", "Greater than or equal to", "Less than or equal to"],
                        key=f"dropdown_{parameter}"
                    )
                selected_operator = comparison_options[dropdown]
                sign = selected_operator
              
        
                    
            
                with col2:
                    limit = st.number_input(
                        f"limit ({parameter})",
                        key=f"lower_{parameter}",
                        value=0.0,
                        step= 0.1

                    )
                with col3:
                    weight = st.number_input(
                        f"Weight ({parameter})",
                        key=f"weight_{parameter}",
                        min_value=0.0,
                        max_value=1.0,
                        step=0.1
                    )
                    
                total_weight += weight
                weights.append(weight)
            
               
            elif types == "Limits":
                with col1:
                    lower = st.number_input(f"Lower limit ({parameter})", key=f"lower_{parameter}", value=0.0)
                with col2:
                    upper = st.number_input(f"Upper limit ({parameter})", key=f"upper_{parameter}", value=1.0)
                with col3:
                    weight = st.number_input(f"Weight ({parameter})", key=f"weight_{parameter}", min_value=0.0, max_value=1.0, step=0.01)
                total_weight += weight
            input_data[parameter] = { types.lower(): { 'lower' : lower if types == "Limits" else limit, 'upper': upper if types == "Limits" else sign, 'weight': weight } }
            st.write("")
            st.write("")

        

        button = st.button( 
            f"Submit {len(st.session_state.get('selected_parameters', [])) } Parameters",
            key=f"submit_{parameter}",
            help=f"Click to submit {parameter} input"
        )

        
        if button:
            if total_weight != 1:
                st.error(f"⚠️ Total weight sum is not equal to 1. Currently: {round(total_weight, 2)}")

           

            elif any(weight <= 0 for weight in weights):
                st.error("⚠️ No weight can be less than or equal to 0. Please adjust the weights.")

            else:
                with st.spinner("Processing..."):
                    if  check(list(input_data.values())):
                        st.session_state["input_data"] = input_data
                        # Simulate processing time
                        time.sleep(2)
                        st.success("✅ Parameters submitted successfully!")
                        st.session_state["home"] = 4
                        st.rerun()

                    else:
                       st.error("⚠️ Invalid input. Please check your parameter values.")







    # Step 5: After submission, show summary
    elif st.session_state.get("home") == 4:
        st.title("Sample File Download")
        st.write("Click the button below to download the sameple on how the columns should be namead.")
        parameters = st.session_state.get("selected_parameters", [])
        input_data = st.session_state.get("input_data", {})
        
        


     


        if file_formation(list(input_data.keys())):
            with open("input_data.csv", "rb") as file:
                st.download_button(
                    label="Download Sample File",
                    data=file,
                    file_name="input_data.csv"
                )


        st.write('')
        st.write('')
        st.write('')
        st.write('')
        st.write('')


        # upload your for the analysis
        uploaded_file = st.file_uploader("Upload your file", type=["csv", "xlsx"], key="file_uploader", help="Upload a CSV or Excel file containing your data for analysis")

        if uploaded_file is not None:
            
            if uploaded_file.name.endswith('.xlsx'):
                file_data = pd.read_excel(uploaded_file)
            else:
                file_data = pd.read_csv(uploaded_file)

            # Check if the uploaded file has the expected columns
            expected_columns = [col.lower() for col in input_data.keys()]
            expected_columns.insert(0, 'id')  # Ensure 'ID' is the first column
            
            
            if check_columns(file_data, expected_columns):
               
                if  'nan' not in [str(name).lower() for name in file_data['ID'].to_list()] :
                    file_name = 'uploaded_file.csv'
                    file_data.to_csv(file_name, index=False)
                    
                    button = st.button(
                        "Submit File",
                        key="submit_file",
                        help="Click to submit your uploaded file for analysis"
                    )
                    if button:
                        with st.spinner("Processing..."):
                            # Simulate processing time
                            time.sleep(2)
                            
                            analysis()
                            st.success("✅ File submitted successfully!")
                            st.session_state["home"] = 5
                            st.rerun()

                else:
                    st.error("⚠️ The 'ID' column in the uploaded file is empty. Please ensure the 'ID' column contains valid data.")
                
            else:
                st.error("⚠️ File is missing some required columns.\n Please Download the sample file and ensure your file has the correct columns.")
                
    elif st.session_state.get("home") == 5:
        st.title("Analysis Results")
        st.write("Click the button below to download the analysis results containing the MPO score data and input data.")
        with open("analysis_results.xlsx", "rb") as file:
            st.download_button(
                label="Download Analysis Results",
                data=file,
                file_name="analysis_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key="download_results"
            )
       
        




elif selected == "Lipinski":
    if "lipinski" not in st.session_state:
        st.session_state["lipinski"] = 1
    if st.session_state.get("lipinski") == 1:
        st.title("Lipinski's Rule of Five")
        st.write(
            "Lipinski's Rule of Five is a set of guidelines to evaluate the drug-likeness of a compound. "
            "It suggests that a compound is more likely to be orally active if it meets the following criteria:"
        )
        st.write('Subtitle: Lipinski\'s Rule of Five')
        st.write(
            """
            1. Molecular weight less than 500 Da
            2. LogP (octanol-water partition coefficient) less than 5
            3. No more than 5 hydrogen bond donors
            4. No more than 10 hydrogen bond acceptors
            5. No more than 1 violation of the above rules is assigned as pass else it fails
            """)
        # input smile as text input or submit a file
        st.write("You can input a SMILES string or upload a file containing SMILES strings for analysis.")
        method = st.radio(
            "Select Input Method",
            options=["SMILES Input", "File Upload"],
            horizontal=True,
            key="input_method"
        )
        if method == "SMILES Input":
            st.write('Input your 1 SMILE   \n  Note: Don\'t enter more than 1 SMILE and You can\'t download the result, you can only View   ')
            smiles_input = st.text_input("Enter SMILES string", key="smiles_input")
        elif method == "File Upload":
            st.write('Ensure you write each SMILE on a newline to get a smooth result \n Note: You can View and Download result')
            file_upload = st.file_uploader("Upload SMILES file", type=["csv", "xlsx"], key="file_upload") 
        
        if st.button("Analyze", key="analyze_button"):
            with st.spinner():
                time.sleep(2)
                if method == "SMILES Input":
                    
                    if len(smiles_input.split()) > 1:
                        st.error(' It takes just one input')

                    if len(smiles_input.split()) == 1:

                        if smiles_input :
                            # Process the SMILES input
                        
                            mol = Chem.MolFromSmiles(smiles_input)
                            if mol is None:
                                st.error("Invalid SMILES string. Please check your input.")

                            else:
                                # Perform Lipinski's analysis
                                MW, HBA, HBD, LogP, scores, Lipinski1 = Lipinski(smiles_input)
                                st.session_state["lipinski"] = 2
                                st.session_state['method'] = method
                                st.rerun()


                                    
                            
                        else:
                            st.error("Please enter a valid SMILES string.")
                elif method == "File Upload":
                    if file_upload is not None:
                        # Read the uploaded file
                        try:
                            if file_upload.name.endswith('.xlsx'):
                                df = pd.read_excel(file_upload)
                            else:
                                df = pd.read_csv(file_upload)
                            df.columns = [col.lower() for col in df.columns]  # Ensure columns are lowercase
                            if 'smiles' in df.columns:
                                smiles_list = df['smiles'].dropna().tolist()
                                ids = df['id'].dropna().tolist()
                                if len(smiles_list) == len(ids):
                                    smiles_input = [smiles for smiles in smiles_list if smiles]
                                else:
                                    st.error("The 'smiles' and 'id' columns must have the same number of entries.")
                            else:
                                st.error("The uploaded file must contain a column named 'smiles' and 'id")
                        
                            if  smiles_input:
                                result = Lipinski(smiles_input)
                                result.insert(0, 'ID', ids) 

                                result.to_csv('lipinski_results.csv', index=False)
                                st.session_state["lipinski"] = 2
                                st.session_state['method'] = method
                                st.rerun()
                            
                        except Exception as e:
                            st.error(f"Error reading file: {e}")
                    else:
                        st.error("Please upload a valid file.")
                    
    if st.session_state.get("lipinski") == 2:
        st.title("Lipinski's Analysis Results")
        if st.session_state.get('method') == 'SMILES Input':
            smiles_input = st.session_state.get("smiles_input", "")
            MW, HBA, HBD, LogP, scores, Lipinski1 = lipi(smiles_input)
            df = pd.DataFrame({
                'SMILES': [smiles_input],
                'Molecular Weight': [MW],
                'HBA': [HBA],
                'HBD': [HBD],
                'LogP': [LogP],
                "Failed Criteria": [scores],
                'Comment': [Lipinski1]
            })
            st.dataframe(df)
        elif st.session_state.get('method') == 'File Upload':
            # Display results for multiple SMILES
            df = pd.read_csv('lipinski_results.csv')
            st.write("Analysis Results:")
            st.dataframe(df)
            st.download_button(
                label="Download Analysis Results",
                data=open('lipinski_results.csv', 'rb').read(),
                file_name='lipinski_results.csv',
                mime='text/csv'
            )


elif selected == "Models":
    st.title("Models")
    st.write("This section is under development")
    st.write("Select a model below to view its details and usage.")   

    st.write("Available Models:")
    radio = st.radio(
        "Select a model",
        options=[None, "PTP1B"],
        key="model_selection",
        horizontal=True
    )

    if radio == "PTP1B":
        st.subheader("PTP1B Model")
        st.write("About the PTP1B model:")
        st.write(
            """
            The Protein Tyrosine Phosphatase 1B (PTP1B) model is designed to predict the inhibitory activity of compounds against the PTP1B enzyme a critical target in diabetes and obesity research.
            It uses a machine learning approach to analyze molecular features and predict the binding affinity of potential inhibitors.
            """

        )

        file = st.file_uploader(
            "Upload a CSV file containing IDS and SMILES for PTP1B analysis and make sure the column is named 'SMILES'",
            type=["csv", "xlsx"],
            key="ptp1b_file_uploader"
        )
        if file is not None:
            if file.name.endswith('.xlsx'):
                df = pd.read_excel(file)
            else:
                df = pd.read_csv(file)

            if 'smiles' in df.columns.str.lower() and 'id' in df.columns.str.lower():
                df.columns = [col.lower() for col in df.columns]
                df.rename(columns={'smiles': 'canonical_smiles', 'id': 'ID'}, inplace=True)
                model = 'models/rf_model.pkl'

                df_new = descriptors(df, model)
                df_new.to_csv('ptp1b_results.csv', index=False)
                # Here you can add the logic to process the SMILES strings with the PTP1B model
            else:
                st.error("The uploaded file must contain 2 columns  named 'ID' and 'SMILES' . Please check your file.")
        else:
            st.warning("Please upload a CSV or Excel file containing SMILES and IDs for PTP1B prediction.")












                        
        

