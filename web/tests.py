from django.test import TestCase
from django.contrib.auth.models import User
from rest_framework.test import APIRequestFactory
from rest_framework.test import force_authenticate

from reads.models import Flowcell, Run, FastqRead, FastqReadType, Barcode
from web.utils import parse_md_cg_pafline


class ParsePafFile(TestCase):

    fixtures = ['fixtures/auxiliary_data.json', ]

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
