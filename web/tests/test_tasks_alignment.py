from django.contrib.auth.models import User
from django.test import TestCase

from alignment.tasks_alignment import calculate_exepected_benefit_2dot0
from reads.models import Flowcell, Run, FastqRead, FastqReadType, Barcode, JobType, JobMaster
from reference.models import ReferenceInfo
from web.utils import parse_md_cg_pafline, parse_mdpaf_alex


class ParsePafFile(TestCase):

    fixtures = ['fixtures/auxiliary_data.json', 'fixtures/reference.json', 'fixtures/yeast.json']

    def test_split_paf_line(self):

        line = "3d8564a8-653e-4dbe-b301-3a2cea209bf8	1129	85	1126	-	ref|NC_001140|	562643	213992	215072	1008	1089	60	NM:i:81	ms:i:1646	AS:i:1646	nn:i:0	tp:A:P	cm:i:76	s1:i:648	s2:i:472	de:f:0.0597	cg:Z:8M1D52M1D55M1D15M4D54M2I28M1I30M1I14M3D23M1D1M1D4M2D53M1D27M2D1M1D47M1D17M1D9M1I25M1D5M1D6M1D66M2D17M1I4M1D4M2D43M1D51M1D3M3D6M1D38M1I22M1I70M1I10M2D14M1D2M1D10M2D26M2D18M1D74M1D25M3D39M1D16M	MD:Z:8^T52^T50C4^A15^CTTG126^CTT23^C1^G4^TT53^T20C2T3^TC1^T23T23^A17^C31T2^C5^A6^A38T0C26^TC18T0T1^C4^GG43^T51^T3^CTT6^G15G0C74T0A0T46^AA1G0A11^C2^T10^TC26^TT18^G0C55T17^A4T2C0T14A1^AAG2T36^A1T14"

        line_tuple = parse_md_cg_pafline(line)

        self.assertEqual(line_tuple.read_id, '3d8564a8-653e-4dbe-b301-3a2cea209bf8')
        self.assertEqual(line_tuple.chromosome, 'NC_001140')

    def test_advanced_mapping(self):

        user = User.objects.create_user(
            username='admin',
            email='admin@minotour.org',
            password='top_secret'
        )

        flowcell = Flowcell.objects.create(
            name='Flowcell 0',
            owner=user
        )
        flowcell.save()

        run = Run.objects.create(
            name='Run 0',
            flowcell=flowcell,
            owner=user
        )
        run.save()

        barcode = Barcode.objects.filter(run=run, name='All reads')[0]

        fastq_read_type_template = FastqReadType.objects.get(pk=1)

        read1 = FastqRead.objects.create(
            run=run,
            flowcell=flowcell,
            read_id='8ad48b94ad4dfc75bdaf3c73485cd99f0dc71eb3',
            read='49869',
            channel='810',
            barcode=barcode,
            barcode_name=barcode.name,
            sequence_length=25,
            quality_average=1.0,
            is_pass=True,
            type=fastq_read_type_template,
            start_time='2018-08-04T11:55:25Z',
            sequence='TCATTTCGTTCGGACATTTATCAGTAGAGCATTAAAAGGTGGAAGATATTGACTCAAAATATTGTGGCATCGCCACGGACAATGGTTCCAGACTATTAAGTACTGATCGACGAGCAGTGTCTGGAGCATTGGTAACTACGGCGATGTTCTTTACTCTGACTACGGTAGTTTATGATGACATAAATCCAGAACAACAGGCATATCGCCAGTATTACTATTGGTATTAGATGACCAGGACAACAGTTGGTTCAATTACCAAAGTGGATAAAAGGATAATGAAAAACTGACCGGGAACTAAATTCCAACTATGCATCGGAAGGTGATAAAGCCGGATGGTGATGGATGGCAATCTGGAGTATAAAATGGTTGAGACGGAAGCGCCAACAGGTTGCACAATTAGTGATGAATATAAGACGATGCAAAAATCACGGATTTACCGCCATTTTGGTCCTAAACGTTACAATTGAAAACACCAGAGCAAATAGTGGCGCGGTTCTTTTAAAGAAGATGGTGTAACGAAGACGCTATCGCAGGGCAAATTCGAATTCAAAATGCAGACGGTACAAAAGTCGCACGAAAACTTAGTATCTAATGCAGATGGAAGCGTCAGAGTAGCAATAGCGCCAGGCGATTGCAATTTGTGAAACAAAGCGTAACTGGCTGTATATTCGTCACAGCGCCAAATTCACGATTGAATTTATTCAAAAATAACTGTATTCACAGCGAAAACTGGTGGTGGTTCCTGAACGAAAGAGATGGAGTGGCTTTCGTGTGGAATTATTCAAAACCGGAACAAAAAGTGTCAAACTTGGTGGACGCAGATGGTAGGTAAACAAACGATTTTAACTGGGAGACTACCAGTTTGCGGCTGGATGCTTCGTTGATTAGTTGTAAACAAGACTGATACGTAACCAAAGAAACTGTGGAATCTGGCGGTTCTGTAAAACAGATGATTGCAAGCTGCGTTGTCTGATGCGAATTTGTTCTGGAACTGCCACTGGAACTAAAGGTAAGGATTTAACAGATGCAAGTGGTGAGAGGATAAGTGGCTGATTTGGCACCGGGCGATTACAAGTTCATTAGAACCAAAGCACTGTGAGCTGCAATTGGATGCCACGCGGACACTTTTTACCGTCAGATTTGGCAATCGATGATCCAAGTGAACAAAGAAGTATTGAACTGGTAGTGTGGTCTTATGACTGAATGACAAATCAAAATTAACAGGTTCAGTTACAAACAAAACGGTGTGGTTTGAAAGATGAAGTAGTGAAGCGATGGTCGACTGCAGGTTGATGGCAGCGCCGGGCGACTGCAATTGAACAGAGCAGAACTCAACAGGTTGTGAGCAGTCTGGCACAGTTAATTCACCATTGAATTTAACCAAAGATGCAGTTCAAGTATTTCAAAACAAATAAAATGTCTTACTGGATCTGTTGTGTTAACTAAGGACTGATGAGCAAAAACTTGCAGACGTTTAAATTGATGAGCAGATAATAACGTCACGAAGAGTCTAACAGCTGGAATAATAATCATGTTGATTGTCGATTAATTGAGGGTTATGAACAATGCCGTTTCGGTCGATGTAATCAGTTTTAATCAAAATCAGTGGCTGGTAACCAAAACAAATATTAAACAGTTCAATGGAATATTCAATTTGTTGATACAAAGAATGTGTGGAATAAGTTCATGCAGGTTTAGTCAGGAAAATATGTAACTAAAGCCAGATCTATCGCGGGTTATAAACTAACAGAAGACCTACAGAATCGGCATACAAAGCAGATCAAAAAGTAACCTTTAACATATGAAAAATAAATCACCATTGTCTGACCCTACTAAATAATTGAAACCTTCCTATTATTAAATACAACTTCAAAACAAGCCACAGATCTTCCTTCAACGAGGAGATAATTCCTAGCTAATCTAATCGTTACGGATTATTCGCTAGTACTAGGACTTTTCTTATAAGAAAAGTCAAAAAGTAAATGATTAATAATTGGTCTTCTGCTGGATAGCGACTTTTTATCATTATATGAAAAAATTCCATATAATGCAGGAAAAGAGGGTTTCATATTGAAAAATGTACATGGTATGCTTTTGTTCACTAAAAATAATCCGCTTTTATTCTTCAAGAACATTTCGACTTCTCGACATTGATATATTATTCGGTCCGTCTTATAAACGAACACATAGGAGAAAACAGCTTAATAAATATATGGTGATGGAGTTTTTATTCGTGGTTACTCGTTTTCAAAACAATATCTTCCATATAATATGACATTGTTGATGAAAATGTAATAATGTAAAATAATCAAAGAGATAAAAGATAGAAAATCAGTAATTAACTAACATGAAATTGAAACAAAACTGGAGAAACAAAGAAGTTTCTATAACTGAACACAGGACAATTATAGAAAATTCACTTTCAGAAGAATATATTGTTAAGAGGCTCCAGGATATATTTAGACAGATGAACGTCACTGGCTGATAAGAGAACTTCTAGTTCAACAAAGAGGCAGAAACCATCATTACCAGAACAGCCTCTAAAAGGGATCTGAAGACAGTTATGTTTCAGATAATATTTTCAAGGTTAAAATGGAACTGGGAGTAGACTTGGTGCATTATGTAAAAATGGCAGACAGTAGTTCCATGAAAATGAATTTTACTTTCTCAGAACTTCTTTAAAGCTGGAGAACAGTTTACCAGATTCATTCAACTTTAGACGGCTACAGACGGAGATTCTTACCTTCAACAATACATGGCA',
            quality="'$$&'%%((,+*/%$%(),)+)+$)%&(&'&%%),+)%&&$'&))),&.2,,++,),,*+(2(%&)'++)*$&%$&()-,0+'%$&&*%('%%%$%(%$$%$&$$%&&%%()%'(&'+,2*$((%(()+.'**,)(%%&''&&%$$&'%%'&%''('')+,'/,1.$'*+%$%%%&(*+(())*+.&('&)))-/+%&'$&&()$%&&%&&'(%&)),+&&&%&&%&(('&%&%.'*%&''''(1''&(%&&)(%*&**)%,(.-*)*(&*()'(&'%')00+($$%*(*,)&*+,/.../+%($$'%'%(&$%&()'*&'&((&'&$%$)&'$$$$$()+/)('%%%*-.)+)'$+)+*0.'&'*+/*)0//%*'(%(%'&'&().*2)%'')+',0((**,000766530/)(((*%%*-.1*((%&(*&,.*(,.,/5-31/*('&%%%*.3-)+,*-..0/,-12/-+.(*(&'),.)'&%'%&&&('*04+076..1(&,+*32(+.(*(-('*''+*)*--10.+%*'&((*.*1.%(&,)%(*..*+/),,+'&*%&''+('(($%#%%(+,*&%%%$++,.'(+*%('),+++'&(*)()%%'$%%&''%$%%&&'''*)(('+'%))+//)%&&(0.(.0'&('&(.,.,,('&'&-'+,%'&$$&&%*%&&'(0+3*%%&%+(&)(-.'%$%$*,+)%$&&$%$'$)%$%%%&%%*,,''&&%'(''(-/.-)('%%&().(&(&(/+'%&&&%$%%$&$%%$%&$$+'$'%$'(('''*-)-+*1.-+',*((,(%&(&%$$%,-.--')''&%&%%&&'.+*+.++''*32-)(&$$'*(**'(%%$%%%%$$#'$$&$$$$$$$%$$(&&%'($%)'$%'(0-)'%$%'&%&%''+*++'.,*%*)&'$$'&&()'')&('&'-5+.%+)',.)*(*,*'+%$$&%&')//.,*%(''(%&$$$)&0/*&()')%%%$%'%''%&&'*')*,10*&&()(&&%''%(&%(')*'%'&%((+*(&&&)*/.'%$&&'()-+16-(/,--3'..,,&,-&++,-/5,&(&%&'&)'--(&'$'$$&'%(''%)'(%+(),,+/0'+((&($$$&(*+*&$'((..'(&.(&''*-1/1'&%%%+,&-*&)&('*.-)%''%'%'%%&'+++,+&'()..0.')'&%%%))'&)(#'%&'&'))+)&%$$$$$()+00%&+)$%%''&*,)&%&*))*%)(0+-)+('**%('&)$&%%%)%&''+.0/..)**'+(()(1'&&%%(')+)+('&%+(()*%)&%)&$%'%&&''(('&&%'-(*)%)$/)$(&'&(%$&%#%(++.*+.3.1*4.,,,,03*,+21')&&$+('%'&%'''-++-***&%&%$/.0/..5.)*+0+('%$'((-,,/*.01-/($'$+'*(,-+'('(,))&%%%$('()&%%&&%'%%%%%*+(%&&)'%%(&&$'%'((&)($$%%''&(&&)'&&'')&&#%$$%)''*'*&')),'%%$%''$-(#&%&&+./2-%$%$$(+1%('('&%'%%&&%(-+)&),+*)*'%''%%%'+13+01,)(,,'+()%'&%%&%%'&(%&'*-+*'(+*)&')*''&).,%%$#%'&&%$%&$'&*0('-,))*%&)(''&('%(&$%((%%'''*0+,0+-'(,-44).00'.**)+,/.1+()'))*)($(%&%++./'-+&*(+-.,'''('&&((*(&**')+/.(&%%%$(*+-)*+)+.1,()&)('()344+)*+-,)((&'&)**((((&&+,-.4*)*%&%%'*(+,%)*261*)()../)..*-%(+'*)+)*&$%%$%$*(*'&%&%%%&'''--('-*%'$&$&$&()+.*)*-,+-,&%+)++.0162,1*&(/0012.-%'-00(/*%((-(*-(+.($%(%&(&''(02//.544+--)*,/(+11.&*-)/.,)$($&%'())0).+,*-(%$&&'((*%$$&&&%%%%&*(*/.-*()(%)&&(-*)*14/.,'*(+))(('&&$$%(.+10,'*)-*-30&$)%+)(,.+,&'$&$%&*')*)($&')++*-40--+,.''&%'-*-(%&)*+($))&*-%%%&&%%)+-*'&'&&)&('%))'('%''&'&'*%&&'(**,,)**+(%&&)'&%(&%$$%$$*))&,01.(''&'$)&'*//+,,'''('')%&(&&*/,)&%%$$$%&%((%%''%)+-,')*(&&'$%$'&&&%&&*(*(''%&%'('','++))(%'((('$$''&,.-(('+.-'++*+,/)'(),,,**'%'%'$),3)-.,&*%$(%&$$&'&%%&$')&+&%))('+/.(-/.-)')*'./3,&)((-14(-&/.1/200**%'$*(*%*)-0/00--3-1*5*()*,2.,&+('(+)01,03-0,.'&'%$$'+,%13++''%')(%&&%&%)&#%$$%)*&)+)../-.-**%&%$$$&&&'&+-2,(*).0)('*(('%'&%%&%&&$&&%'&'(+-))$*(+)&*+.((*())(%%&'%'&'&.(%&$'('$$(('(/**,(''%%)&%&&(+.(('*''%&')%&%'%%'$'%(&'&''*'%&%&&++*+(&%$%&(''%(%,.,)$&&$(*,**()(.44.('+,+)&&%%(&%%.%&%$')('),%%&#%**',0,'$%$&&),/,,,3,**().(%'''%(((%%&&%()+''%.((('%)(*+%%&$%+()*%'%%&%",
        )
        read1.save()

        print(Flowcell.objects.all())

    def test_advanced_mapping2(self):

        JOB_TYPE_READUNTIL=15
        REFERENCE_INFO_S_CER_S288C=18

        flowcell = Flowcell.objects.get(pk=1)

        run = Run.objects.get(pk=1)

        job_type = JobType.objects.get(pk=JOB_TYPE_READUNTIL)

        reference_info = ReferenceInfo.objects.get(pk=REFERENCE_INFO_S_CER_S288C)

        job_master = JobMaster.objects.create(
            run=run,
            flowcell=flowcell,
            job_type=job_type,
            reference=reference_info
        )

        calculate_exepected_benefit_2dot0(flowcell.id, job_master.id)

    def test_parse_mdpaf_alex(self):

        line = "0	tp:A:P	cm:i:652	s1:i:4465	s2:i:0	de:f:0.0441	rl:i:0	cg:Z:6M1I9M2I8M1D17M1D2M1D12M1I29M1D42M1I6M1I29M1D18M1I11M2I6M1I3M1D33M1I10M1I96M1I6M2I155M1D38M1I2M2I14M1D35M2D32M1D79M1D12M3D60M1I5M1I65M1I88M1I6M1I23M1I37M2D26M1I10M1I91M1D26M2D24M1I13M1I6M2I29M1D11M1I9M1I83M1I4M4D5M1I4M1I44M1D18M1I16M1D13M1D82M1I8M1D10M1D29M2I10M2D45M1D6M1D11M1I7M2D34M1I23M2D42M1D18M1D48M2D9M1D4M8D27M2D121M3D25M1D25M1D26M1I52M1D10M1I9M1I11M1D4M3D27M1D42M1I11M1D8M1I24M1D69M1D31M2D9M2D19M2D99M1D22M1D11M1I101M1I14M1D14M2D40M3D43M1D1M1D47M3D1M2D29M1D54M1D12M1D5M1D1M1D20M1I55M1D28M1D4M1D58M1D15M1D80M1D36M1D24M1D2M1I47M1I39M3D76M1I39M1D62M2I4M1I48M2D2M1D14M2I19M1D1M1D26M2D10M1I64M1D24M3D34M1D142M4D109M1D9M2D29M3D4M1I46M1I2M1I51M1D11M1I71M1D4M1I19M1I70M1D33M1D19M1I22M1I53M1D44M1D2M1D11M1D26M1I38M2D26M1D11M1D43M2D21M1D31M2D52M3D14M1D10M2D123M1D80M2D15M1I44M1D53M1D9M1I14M3D21M1D12M1I7M1D65M1D13M1D173M1D97M1D40M1I135M1I22M1D6M1D31M1I16M1I61M1I33M1I70M1I12M1D27M1D18M1I22M3I85M1I35M	MD:Z:6G14A1^T17^A2^G41^A75A1^G31G1A4^C112T0G0C0T0T183^T36T5C11^A35^AA32^A79^A12^AAT219G0A60G0G1^TT13A113^C26^AC14A22T34^T107^CAAT53^A1C17G14^T13^A19A1T2G57A0T6^C10^G30C5G2^GC45^A0T5^C18^CT57^CC42^C18^A0A47^AG9^C4^GAAATGAC1C25^GA121^CGA25^C25^G31A1T44^A30^T4^TTT27^T53^T13A0A17^A69^A31^TC9^GA19^CA10G1A0T0C66G0A16^T22^C1G43G80^T14^CC40^GTT43^C1^T44A2^CAA1^CT29^G0C53^T12^C3A1^T1^C75^T28^C4^C0A57^T12G2^G24A55^A36^A24^A88^GTT115^A114^AA0A1^G2G30^T1^G26^GA10G25A21A15^A24^TGA34^G0C64G0A75^TTGG109^T0T1T6^AT29^GGT1A2A72A2G22^G4A0T76^C38G54^A33^A94^T44^A2^A0A10^A9G16G37^GA26^T11^A43^AT21^A31^TA0T51^TGC14^C10^GT11T46C0T63^T80^TA15G43^T53^C1T7G13^GGC19G1^C17G1^A39A0G24^G0C12^A18C31G122^A97^T38C6G101T0T1G0A45^T6^T223^T27^C160"

        mismatcharray, matcharray, mapstart, mapend, maporientation, reference, referencelength, readlength = parse_mdpaf_alex(
            line)

        print(mismatcharray)
        print(mismatcharray.dtype)
        print(matcharray)
        print(matcharray.dtype)

