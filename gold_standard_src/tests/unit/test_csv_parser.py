from gold_standard_src.gold_standard.parsers.csv_parser import csv_corvar_to_num


def test_csv_corvar_to_num():
    csv_corvar = "------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "-------------------------------------k----------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "-------------------------------------em-----p---------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "-------------------------q------------------------------p---" \
        "-----------------------k------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "----------------------------------------------t-----------f-" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "-------------------------------------------g----------------" \
        "-------------------------e-------------------lk-------------" \
        "------------n---------------------------------------l-p-----" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "-----------------------------------------l-----l------------" \
        "------------------------------------------------------------" \
        "----------------n-------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "---------------------------------------------------------t--" \
        "-----------------------------------------------------------d" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "-------------kp---------------------------------------------" \
        "---------------------------v-------------------------q------" \
        "------------------------------------------------------------" \
        "----a-------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "----------------------------------------------l-------------" \
        "--------------------------m-k-------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "--------------------------------------------------ia--------" \
        "------------------------------------------------------------" \
        "-------d-----------------------------------------e----------" \
        "-l----------------------------------------g-----------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------ --EIFKFEA -" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "pgr-- VTRYLSSQRLIKEACDES- ---------- ------ ----------------" \
        "------------------------------------------------------------" \
        "------- ---------- ------------------------r---f------dk---n" \
        "--l--sq-------------a-l------k---------------fv-r--d-------f" \
        "a------------gdg--------------lf--------------------t-------" \
        "----------swt----------------h--e------kn----------------- -" \
        "WKKAHNILLPSFSQ ---------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------ QAMKGYHAMMV" \
        "DIAVQLVQKWE- ----------------r----------------l--------n----" \
        "a---------d-------------------------------------------------" \
        "--------------------------- -EHIEVPE- ----------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "--------------------------------------------------------- DM" \
        "TRLTLDTIGLCGFN ---------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------ -----------" \
        "---------- ----------------------------------------------y--" \
        "---------------------r------------f--n---s------f---------y-" \
        "----------r-------------------------d--------------q--------" \
        "---------p-hpfi------------ts------mv-----------------------" \
        "---ralde---a-----------------m-----nk-------------l---------" \
        "------------n---------------p------------------------d------" \
        "--------------------------------------d------------------p-a" \
        "----------------------------y---d--------e-n----------------" \
        "-----------k----------rqfqe---------------di--------kv------" \
        "---------------m--------------n-----d---------------l-----v-" \
        "-d------------------k-----------------------------i---------" \
        "-------i------------ad-----r--------------------------------" \
        "----------k-----------a----s--------------------------------" \
        "---------------------------------------------ge-------------" \
        "-------q--------------s---------------------------- -DDLLTHM" \
        "LNGKD- -------------p---------------e--------t--------------" \
        "------------------------------------------------------------" \
        "------------- GEPLDDENIRYQIITF -----------------------------" \
        "--------------------- LIAGHE -------------------------------" \
        "-- TTSGLLSFALYFLVKN ----------------------------------------" \
        "--------------- PHVLQKAAEEA -arv--------l-----v-------------" \
        "--------------------------------------d--------p------------" \
        "vp--------sy------------k------------------q---vkq----------" \
        "-----l-------------------------- KYVGMVLNEALRLWPTA ---------" \
        "---------------------------------------------p- AFSLYAKEDTVL" \
        "G -g--------------------------------------------------------" \
        "--------- EYPLEKGDELMVLIPQL ------------------------ HRDKTIW" \
        "G -d-------------------- DVEEFRPER-- --f--------------------" \
        "-----en--------------------p--------------------------------" \
        "-----------------------s---------------------------a--------" \
        "--------i---p---------------------qh------------------------" \
        " -AFKPFGN -------------------------------------- GQRACIGQQFA" \
        "LHEATLVLGMMLKHFDFEDHT -----------n------------y---------e---" \
        "------------------------------------------------------------" \
        " LDIKETLTLKPE ---------------------------- GFVVKAK -sk------" \
        "----------------------------------------------------------k-" \
        "-----------------i----------pl------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "---"
    print csv_corvar_to_num(csv_corvar)
