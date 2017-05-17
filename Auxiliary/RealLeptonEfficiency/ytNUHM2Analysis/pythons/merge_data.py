#!/usr/bin/python
import os, ROOT, getopt

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "b", ["skimmed", "RLE", "RLE-trigger="])
    except getopt.GetoptError:
        print 'merge_data.py -b --skimmed/--RLE' # -b: batch mode
        os.sys.exit(2)

    for opt, arg in opts:
        if opt == "--skimmed":
            merge_skim_data_hist()
        elif opt == "--RLE":
            merge_RLE_data_hist("")
        elif opt == "--RLE-trigger":
            trigger = arg
            merge_RLE_data_hist(trigger)
            #merge_RLE_data_hist("dilepton_trigger")
            #merge_RLE_data_hist("dilepton_trigger_but_fail_single_lepton_trigger")

def merge_skim_data_hist():
    command = "hadd hist-merged-Data.root "
    for folder in sorted(os.listdir("./")):
        if os.path.isdir(folder) and "skim_Data_merged" in folder:
            for file in sorted(os.listdir(folder)):
                if file.startswith("hist-") and file.endswith("root"):
                    #print folder + "/" + file
                    command = command + folder + "/" + file + " "
    print command
    os.system(command)

def merge_RLE_data_hist(trigger):
    if trigger:
        trigger = "_" + trigger

    command = "hadd hist-RLE-merged-data-elec" + trigger + ".root "
    for folder in sorted(os.listdir("./")):
        if os.path.isdir(folder) and "RLE_Data_merged_data" in folder and folder.endswith("electron" + trigger):
            for file in sorted(os.listdir(folder)):
                if file.startswith("hist-") and file.endswith("root"):
                    #print folder + "/" + file
                    command = command + folder + "/" + file + " "

    print command
    os.system(command)

    command2 = "hadd hist-RLE-merged-data-muon" + trigger + ".root "
    for folder in sorted(os.listdir("./")):
        if os.path.isdir(folder) and "RLE_Data_merged_data" in folder and folder.endswith("muon" + trigger):
            for file in sorted(os.listdir(folder)):
                if file.startswith("hist-") and file.endswith("root"):
                    #print folder + "/" + file
                    command2 = command2 + folder + "/" + file + " "

    print command2
    os.system(command2)

    #!/usr/bin/python
    #import os, ROOT
    #el_list = []
    #mu_list = []
    #for file in sorted(os.listdir("./")):
    #   if os.path.isdir(file) and "RLE_Data" in file:
    #       if "electron" in file:
    #           for root in os.listdir(file):
    #               if root.startswith("hist-") and root.endswith(".root"):
    #                   el_list.append(os.getcwd() + "/" + file + "/" + root)
    #       if "muon" in file:
    #           for root in os.listdir(file):
    #               if root.startswith("hist-") and root.endswith(".root"):
    #                   mu_list.append(os.getcwd() + "/" + file + "/" + root)
    #
    #hadd_el = "hadd hist_all_data_electron.root "
    #for item in el_list:
    #   hadd_el += item + " "
    #
    #hadd_mu = "hadd hist_all_data_muon.root "
    #for item in mu_list:
    #   hadd_mu += item + " "   
    #
    #os.system(hadd_el)
    #os.system(hadd_mu)

if __name__ == "__main__":
    main(os.sys.argv[1:])
