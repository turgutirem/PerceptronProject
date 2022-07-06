using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

/*
 çok katmanlı yapıda aynı doğrudaki verileri ayrıştımaya çalışmak için çok katmana ihtiyacımız var! tek katmanda linear olarak 
çalışoupruz. 2 sınıf için ilk katmanda kendimşiz belirleyerek nöron olşuturuyoz,  gizli katman için belirlenen nöron sayısı da gizli  katman için belirlenen 
nöron sayısı oluyor!
artık farklı bir uzayda işlem yapıyoruz  yani çözüm uzayında eğrilerle belirtiyoruz.
yani linear olarak ayrıştıramazsak çok katmana geçtik.
gizli katman nöron sayısı kadar doğru
nöron sayısını arttırdığımızda hatayı daha küçültebilirz diye yrumladım.
en az sınıf sayısı kadar nöron koymamız lazımmm
hata azalmıyor aksine büyüyorsa bu yapının bunu çözmekte yetersiz kaldığı anlamına gelir
 */
//sınıf sayısı kadar nöron
//hata değeri ya da max iterasyon (iterasyon sayısı) sonlandırma kriteri olarak kullanılır.
// Binary aktivasyon fonk. ayrık değer alır.
//Bİnary'de tüm değerler doğru sınıflandırılırsa hata 0 olur. Süreklide öyle değil.
//Gizli katman uzayın taşınması gibidir. Yani iki boyutludan üç boyutluya geçiyormuş gibi.
//Gizli katmandaki nöron sayısı doğru sayısına eşittir.
//perceptron öğrenme kuralı -> delta w = c (d-o).x  
// c = öğrenme katsayısı, d = gerçek sınıf etiketi , o = hesaplanan sınıf etiketi, x = input
/*Okulda sınıflandırma problemleri hakkında çalışma yapacağız.
2 sınıf için tek bir nöron yeterlidir.Bunu nedeni zaten sonuca göre ya 1. sınıfa aittir ya da 2 seçeneğimiz bu kadar.
Fakat 2 den fazla sınıflar için o kadar(sınıf sayısı kadar) nöron olması gerekir.
perception öğrenmeye bakıcaz.^w=c(d-f(net)).x formüldür.f(net)->hesaplanan etiket
Yukarıda c=öğrenme katsayısı,d=grçek sınıf etiketi.Burada ağırlık güncellemesi yapıyoruz.
Hata fonksiyonu : e=1/2(d-f(net))^2 hatayı hesaplıyoruz ve hatanın azalmasına bağlı olarak ağı eğitiyoruz.
hesaplamalarımıza bias'ı 1 alarak devam edicez.*/
//w1=w1+c(d-o)*x1..... devaö ediyor. sonunda wn=wn+c(d-o)*xn olacak. en son w0=w0+c(d-o)*bias kolaylık olsun diye w0 sonda gibi işlem yaptık
//sürekli değerleri normalize yapmdan kullanmıyoruz.
namespace PerceptronProject//perceptron projesi
{

    class Data//veri modeli tanımlıyoruz.
    {
        public double[] input { get; set; }
        public int output { get; set; }

        public Data(double[] input, int output)//x ler double veri tipinde çünkü matematiksel işlemler var,ağırlıklar vs.
                                               //çıkış değeri de int çünkü sınıfların etiket karşılığıdır.0,1,2 gibi 
        {
            this.input = input;
            this.output = output;
        }
    }

    class Neuron
    {
        public double bias { get; set; }
        public double[] w { get; set; }

        public Function function { get; set; }
        public Neuron(int dimension, double bias, Function function)
        {
            w = new double[dimension];
            for (int i = 0; i < dimension; i++)
            {
                w[i] = new Random().NextDouble();
            }
            this.bias = bias;

            this.function = function;
        }

        public String getW()
        {
            String data = "";

            for (int i = 0; i < w.Length; i++)
            {
                data = data + "" + w[i];
            }

            return data;

        }
        public double getClass(int classNumber, int classIndex)
        {
            if (classNumber == classIndex)
                return 1;
            else return -1;
        }
    }

    abstract class Function
    {
        public double net(double[] input, double[] w, double bias)
        {
            double sum = 0;

            for (int i = 0; i < input.Length; i++)
            {
                sum += input[i] * w[i];
            }

            return sum + w[w.Length - 1] * bias;
        }

        public abstract double calculate(double net);
    }
    //verileri ayırdığı an sonlanıyor
    class Binary : Function // bu da aktivason fonk çeşidi hata 0 ın altına düştüğü an işlem sonlanır
    {
        public override double calculate(double net)
        {
            if (net > 0)
                return 1;
            return -1;
        }
    }
    //burada belirlenen hatanın en altına düşünceye kadr buluyor normalizeden sonra yani nokataların hepsine eşit doğrular çizmeye çalışıyoruz.
    //en iyi eğer 0,01 belirlenmişse onun altına düştüğü an sonlanır 
    //neden çok katmanlı yapı kllanıyoruz tek nöronda bu şöyle oluyor çizilen doğru kendi verilerini üst taraafta tanımlıyor altında kalan değerler ise başka bir sınıfa aitmiş gibi oluyor yani diğer veriler iç içe yan yana olmasını umursamıyor diyebiliriz.
    //tek katmanda 2 den fazla sınıf olursa doğru sayımız sınıf sayısı kadar olur.
    //her zaman hata var minimize etmeye çalışıyor
    //burda mutlaka normalize edicez
    class Continuous : Function //aktivasyon ya da transfer fonksiyonları burada her zaman bir hata değeri olacak sonlanmaz!
    {
        public override double calculate(double net)
        {
            return 1 / (1 + Math.Pow(Math.E, -net));
        }
        public double derivateCalculate(double net) => 1.0 - Math.Pow(net, 2);



    }

    public partial class Form1 : Form
    {
        private List<Data> dataList;
         int classCount;

        public Form1()
        {
            InitializeComponent();
            dataList = new List<Data>();

            double[] d1 = { 1, 2 };
            double[] d2 = { 2, 3 };
            double[] d3 = { -1, 0 };
            double[] d4 = { -1, -2 };
            double[] d5 = { 1, 0 };
            double[] d6 = { -1, 1 };

            dataList.Add(new Data(d1, 1));
            dataList.Add(new Data(d2, 1));
            dataList.Add(new Data(d3, -1));
            dataList.Add(new Data(d4, -1));
            dataList.Add(new Data(d5, 1));
            dataList.Add(new Data(d6, 1));
            HashSet<int> classArray = new HashSet<int>();
            for (int i = 0; i < dataList.Count; i++)
            {
                classArray.Add(dataList[i].output);
                classCount = classArray.Count;
            }

        }

        private void button1_Click(object sender, EventArgs e)
        {
            singleLayerSingleNeuron(0.1, (dataList[0].input.Length + 1), 1);
            singleLayerMultiNeuron(0.1, (dataList[0].input.Length + 1), 1);
        }

        public void singleLayerSingleNeuron(double c, int dimension, double bias)
        {
            Neuron neuron = new Neuron(dimension, bias, new Binary());

            while (true)
            {
                double error = 0;
                for (int i = 0; i < dataList.Count; i++)
                {
                    double net = neuron.function.net(dataList[i].input, neuron.w, neuron.bias);

                    double fnet = neuron.function.calculate(net);

                    for (int j = 0; j < (dimension - 1); j++)
                        neuron.w[j] = neuron.w[j] + c * (dataList[i].output - fnet) * dataList[i].input[j];

                    neuron.w[dimension - 1] = neuron.w[dimension - 1] + c * (dataList[i].output - fnet) * neuron.bias;

                    error = error + Math.Pow(dataList[i].output - fnet, 2) / 2;
                    

                  
                }
                if (error < 0.1)
                    break;
                MessageBox.Show(neuron.getW());
            }

            //Console.WriteLine(neuron.getW());
        }

        public void singleLayerMultiNeuron(double c, int dimension, double bias)
        {
            
            List<Neuron> neuronList = new List<Neuron>();
            for (int i=0;i<classCount;i++)
            {
                neuronList.Add(new Neuron(dimension, bias, new Binary()));
            }
            while (true)
            {
                double error = 0;//1 itersyon
                for (int i=0;i<dataList.Count;i++)
                {
                    for (int j=0;j<classCount;j++)
                    {
                        Neuron neuron = neuronList[j];
                        double net = neuronList[j].function.net(dataList[i].input, neuron.w, neuron.bias);//net hesabı yapıldı
                        double fnet = neuron.function.calculate(net);//fnet hesabı yaptık

                        for (int k = 0; k < (dimension - 1); k++)
                            neuron.w[k] = neuron.w[k] + c * (neuron.getClass(dataList[i].output,j) - fnet) * dataList[i].input[k];

                        neuron.w[dimension - 1] = neuron.w[dimension - 1]+ c * (neuron.getClass(dataList[i].output,j) - fnet) * neuron.bias;
                        error = error + Math.Pow(neuron.getClass(dataList[i].output,j) - fnet, 2) / 2;//1/2*(d-fnet )karesi
                    }
                }
                if (error < 0.1*classCount)
                    break;
            }
        }
        //o sınıfla eşleiyor mu eşleşmiyor muuuu?
        int maxIteration = 1000;
        public void multiLayerMultiNeuron(double c, int dimension, double bias,int neuronSize) 
        {
            List<Neuron> inputNeuronList = new List<Neuron>();
            List<Neuron> hiddenNeuronList = new List<Neuron>();
            //giriş katmandaki ağırlıklar normal v gizli katman w
            double[] fnetV = new double[neuronSize];//giriş katmanın çıkışı artık fnetv oldu.

            for (int i = 0; i < classCount; i++)//nöronları oluşturduk //sınıf sayısı kadar
            {
                inputNeuronList.Add(new Neuron(dimension, bias, new Continuous()));
            }
            for (int i = 0; i < neuronSize; i++)//nöronları oluşturduk //nöron boyutu kadar
            {
                hiddenNeuronList.Add(new Neuron(dimension, bias, new Continuous()));
            }
            int iteration = 0;
            while (true)//iterasyon döngüsü              //uzay dönüşümünde 2 boyuttan 3 boyutlu bir uzaya geçiyor ve daha çok doğruya ihtiyac duyuyoruz.
            {
                iteration = iteration + 1;
                if (iteration > maxIteration)
                {
                    break;
                }
                for (int i=0;i<dataList.Count;i++)
                {
                    for (int j = 0; j < classCount; j++)//gizli katmanımızdaki nöronları döndüren döngü
                    {
                        Neuron hiddenNeuron = hiddenNeuronList[j];
                        for (int k = 0; k < neuronSize; k++)//giriş katmanındaki nöronlar için tüm işlemler yapılıyor
                        {
                            Neuron inputNeuron = inputNeuronList[k];
                            double net = inputNeuron.function.net(dataList[i].input, inputNeuron.w, inputNeuron.bias);//net hesabı yapıldı
                            double fnet = inputNeuron.function.calculate(net);//fnet hesabı yaptık
                            fnetV[k] = inputNeuron.function.calculate(net);//giriş katmanındaki tüm fnetler hesaplandı ve dizide saklandı
                        }
                        double netW = hiddenNeuron.function.net(fnetV, hiddenNeuron.w, hiddenNeuron.bias);//toplama fonk çalışıyor gizli katmandaki
                        double fnetW = hiddenNeuron.function.calculate(netW);//sonra aktivasyon fonka netw gönderilip fnet hesaplanıyor


                        for (int k = 0; k <neuronSize; k++)//ağırlık güncellemesi gizli katmandki perceptron kuralı ile
                            hiddenNeuron.w[k] = hiddenNeuron.w[k] + c * (hiddenNeuron.getClass(dataList[i].output, j) - fnetW) * fnetV[k];
                        //buradan sonra geri yayılım ile bir önceki katmanın(giriş katman olabilir) ağırlıklarını güncelliyoruz ama tüm nöronların tek bir nöron değil
                        hiddenNeuron.w[neuronSize] = hiddenNeuron.w[neuronSize] + c * (hiddenNeuron.getClass(dataList[i].output, j) - fnetW) * hiddenNeuron.bias;//perceptron kuralı

                        Continuous cont = (Continuous)hiddenNeuron.function;//türev iiçin uğraşıyoruzz//polimoorfizm
                        for (int k = 0; k < neuronSize; k++)//giriş katmanındaki nöronlar için tüm işlemler yapılıyor
                        {
                            Neuron inputNeuron = inputNeuronList[k];
                            for (int m = 0; m < dimension; m++)
                                inputNeuron.w[m] = inputNeuron.w[m] + c * (hiddenNeuron.getClass(dataList[i].output, j) - fnetW)//öğrenme kuralı perceptron
                                    * cont.derivateCalculate(fnetV[k])*dataList[i].input[m];//nöron bilgisinde işlem yapıyorsak onu gönderdik fnetv k da

                            //burada sadece aktivasyon  kuralının türevi ilave edilerekk işlem yapılııyor. d-o işlemi c den sonraki işlem 

                            inputNeuron.w[dimension] = inputNeuron.w[dimension] + c * (hiddenNeuron.getClass(dataList[i].output, j) - fnetW)//öğrenme kuralı perceptron
                                    * cont.derivateCalculate(fnetV[k]) * inputNeuron.bias;//nöron bilgisinde işlem yapıyorsak onu gönderdik fnetv k da
                        }


                      /*  foreach (Neuron inputNeuron in inputNeuronList)
                        {


                            for (int k = 0; k < dimension; k++)
                            {
                                inputNeuron.w[k] = inputNeuron.w[k] + c * (inputNeuron.getClass(dataList[i].output, j) - fnetW)* cont.derivateCalculate() ;

                            }


                        }*/

                        //double constant = cFuntion.calculate(c,list[i].getClass(),j,fnetW);
                        //cFunction.derivateCalculate(fnetV[k])
                       
                       // cont.derivateCalculate();
                       
                    }
                }
            }
        }
    }
}

