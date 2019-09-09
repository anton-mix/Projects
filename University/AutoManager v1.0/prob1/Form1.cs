using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Xml;

namespace prob1 {
    public partial class Form1 : Form {
        int size = 0;
        //Нормированные абсолютные координаты
        private const int MOUSEEVENTF_ABSOLUTE = 0x8000;
        //Нажатие на левую кнопку мыши
        private const int MOUSEEVENTF_LEFTDOWN = 0x0002;
        //Поднятие левой кнопки мыши
        private const int MOUSEEVENTF_LEFTUP = 0x0004;
        //перемещение указателя мыши
        private const int MOUSEEVENTF_MOVE = 0x0001;
        //Нажатие на левую кнопку мыши
        private const int MOUSEEVENTF_RIGHTDOWN = 0x0008;
        //Поднятие левой кнопки мыши
        private const int MOUSEEVENTF_RIGHTUP = 0x0010;

        [System.Runtime.InteropServices.DllImport("user32.dll", CharSet = System.Runtime.InteropServices.CharSet.Auto,
    CallingConvention = System.Runtime.InteropServices.CallingConvention.StdCall)]
        public static extern void mouse_event (uint dwFlags, int dx, int dy, int dwData, int dwExtraInfo);
        public Form1 () {
            InitializeComponent();
            comboBox1.Items.Add("Кнопки мыши");
            comboBox1.Items.Add("Перемешение мыши");
            comboBox1.Items.Add("Клавиатура");
            comboBox1.Items.Add("Ожидание");
            comboBox1.SelectedIndex = 0;
            toolStripStatusLabel2.Text = "X = " + SystemInformation.PrimaryMonitorSize.Width + " Y = " + SystemInformation.PrimaryMonitorSize.Height;
            LoadName();
            numericUpDown2.Maximum = size;
        }
        private void comboBox1_SelectedIndexChanged (object sender, EventArgs e) {
            comboBox2.Text = "";
            switch (comboBox1.SelectedIndex) {
                case 0:
                    comboBox2.Items.Clear();
                    comboBox2.DropDownStyle = ComboBoxStyle.DropDownList;
                    comboBox2.Items.Add("Правая кнопка мыши");
                    comboBox2.Items.Add("Левая кнопка мыши");
                    comboBox2.Items.Add("Правая кнопка мыши вниз");
                    comboBox2.Items.Add("Правая кнопка мыши вверх");
                    comboBox2.Items.Add("Левая кнопка мыши вниз");
                    comboBox2.Items.Add("Левая кнопка мыши вверх");
                    break;
                case 1:
                    comboBox2.Items.Clear();
                    comboBox2.DropDownStyle = ComboBoxStyle.Simple;
                    break;
                case 2:
                    comboBox2.Items.Clear();
                    comboBox2.DropDownStyle = ComboBoxStyle.Simple;
                    break;
                case 3:
                    comboBox2.Items.Clear();
                    comboBox2.DropDownStyle = ComboBoxStyle.DropDown;
                    for (int i = 1; i <= 12; i++) 
                        comboBox2.Items.Add(i * 5 + " сек");
                    break;
            }
            button2.Enabled = (comboBox1.Text == "" || comboBox2.Text == "") ? false : true;
        }
        private void comboBox2_SelectedIndexChanged (object sender, EventArgs e) {
            button2.Enabled = (comboBox1.Text == "" || comboBox2.Text == "") ? false : true;
        }
        private void button2_Click (object sender, EventArgs e) {
            if (comboBox1.Text == "Ожидание" && comboBox2.Text.Split(' ').Length == 1)
                comboBox2.Text = comboBox2.Text + " сек";
            if(numericUpDown2.Value == 0)
                listBox1.Items.Add("Тип команды: " + comboBox1.Text + ", команда: " + comboBox2.Text);
            else {
                string[] buff = new string[listBox1.Items.Count];
                for (int i = 0; i < listBox1.Items.Count; i++) {
                    buff[i] = listBox1.Items[i].ToString();
                }
                listBox1.Items.Clear();
                for (int i = 0; i < numericUpDown2.Value; i++) {
                    listBox1.Items.Add(buff[i]);
                }
                listBox1.Items.Add("Тип команды: " + comboBox1.Text + ", команда: " + comboBox2.Text);
                for (int i = (int)numericUpDown2.Value; i < buff.Length; i++) {
                    listBox1.Items.Add(buff[i]);
                }
                listBox1.SelectedIndex = (int)numericUpDown2.Value - 1;
            }
            size++;
            numericUpDown2.Maximum = size;
            button3.Enabled = button1.Enabled = true;
            
        }
        private void button3_Click (object sender, EventArgs e) {
            if (listBox1.SelectedIndex != -1) {
                listBox1.Items.RemoveAt(listBox1.SelectedIndex);
                size--;
                numericUpDown2.Maximum = size;
                if (size <= 0) {
                    size = 0;
                    button3.Enabled = button1.Enabled = false;
                }
            }
        }
        private void comboBox2_TextChanged (object sender, EventArgs e) {
            switch (comboBox1.SelectedIndex) {
                case 1: {
                        string[] split = comboBox2.Text.Split(' ');
                        int x;
                        if (split.Length == 2 && int.TryParse(split[0], out x) && int.TryParse(split[1], out x) ) {
                            split[0] = (Convert.ToInt32(split[0]) <= SystemInformation.PrimaryMonitorSize.Width) ? split[0] : SystemInformation.PrimaryMonitorSize.Width.ToString();
                            split[1] = (Convert.ToInt32(split[1]) <= SystemInformation.PrimaryMonitorSize.Height) ? split[1] : SystemInformation.PrimaryMonitorSize.Height.ToString();
                            comboBox2.Text = split[0] + " " + split[1];
                            comboBox2.SelectionStart = comboBox2.Text.Length;
                            button2.Enabled = true;
                        }
                        else
                            button2.Enabled = false;
                    }
                    break;
                case 2:
                    button2.Enabled = (comboBox2.Text.Length > 0)? true : false;
                    break;
                case 3: {
                        bool prov = false;
                        int x;
                        string[] split = comboBox2.Text.Split(' ');
                        switch (split.Length) {
                            case 1:
                                prov = int.TryParse(split[0], out x);
                                break;
                            case 2:
                                if (int.TryParse(split[0], out x) && split[1] == "сек")
                                    prov = true;
                                break;
                        }
                        button2.Enabled = prov;
                    }
                    break;
            }
        }
        private void timer1_Tick (object sender, EventArgs e) {
            toolStripStatusLabel1.Text = "X: " + Cursor.Position.X + " Y: " + Cursor.Position.Y;
        }
        private void button1_Click (object sender, EventArgs e) {
            for (int z = 0; z < numericUpDown1.Value; z++) 
                for (int i = 0; i < listBox1.Items.Count; i++) {
                    listBox1.SelectedIndex = i;
                    //numericUpDown2.Value = i;
                    //numericUpDown2.Invalidate();
                    listBox1.Invalidate();
                    string[] split = listBox1.SelectedItem.ToString().Split(' ');
                    switch (split[2] + " " + split[3]) {
                        case "Кнопки мыши,":
                            switch (split[5] + " " + split[6] + " " + split[7]) {
                                case "Правая кнопка мыши":
                                    mouse_event(MOUSEEVENTF_ABSOLUTE | MOUSEEVENTF_RIGHTDOWN | MOUSEEVENTF_RIGHTUP, Cursor.Position.X, Cursor.Position.Y, 0, 0);
                                    break;
                                case "Левая кнопка мыши":
                                    mouse_event(MOUSEEVENTF_ABSOLUTE | MOUSEEVENTF_LEFTDOWN | MOUSEEVENTF_LEFTUP, Cursor.Position.X, Cursor.Position.Y, 0, 0);
                                    break;
                                case "Правая кнопка мыши вниз":
                                    mouse_event(MOUSEEVENTF_ABSOLUTE | MOUSEEVENTF_RIGHTDOWN, Cursor.Position.X, Cursor.Position.Y, 0, 0);
                                    break;
                                case "Правая кнопка мыши вверх":
                                    mouse_event(MOUSEEVENTF_ABSOLUTE | MOUSEEVENTF_RIGHTUP, Cursor.Position.X, Cursor.Position.Y, 0, 0);
                                    break;
                                case "Левая кнопка мыши вниз":
                                    mouse_event(MOUSEEVENTF_ABSOLUTE | MOUSEEVENTF_LEFTDOWN, Cursor.Position.X, Cursor.Position.Y, 0, 0);
                                    break;
                                case "Левая кнопка мыши вверх":
                                    mouse_event(MOUSEEVENTF_ABSOLUTE | MOUSEEVENTF_LEFTUP, Cursor.Position.X, Cursor.Position.Y, 0, 0);
                                    break;
                            }
                            break;
                        case "Перемешение мыши,":
                            Cursor.Position = new Point(Convert.ToInt32(split[5]), Convert.ToInt32(split[6]));
                            break;
                        case "Клавиатура, команда:":
                            string word = "";
                            for (int j = 0; j < split.Length - 4; j++)
                                word = word + split[4 + j] + " ";
                            SendKeys.Send(word);
                            break;
                        case "Ожидание, команда:":
                            long tiks = DateTime.Now.Ticks + Convert.ToInt32(split[4]) * 10000000;
                            while (true)
                                if (DateTime.Now.Ticks >= tiks)
                                    break;
                            break;
                    }
                }
            
            listBox1.SelectedIndex = -1;
        }

        private void comboBox2_KeyUp (object sender, KeyEventArgs e) {
            if (e.KeyCode == Keys.Enter && button2.Enabled)
                button2.PerformClick();
        }
        bool LoadName () {
            XmlDocument doc = new XmlDocument();
            try {
                doc.Load(Environment.CurrentDirectory + "\\Save.xml");
            }
            catch {
                return false;
            }
            comboBox3.Items.Clear();
            for (int i = 0; i < doc.ChildNodes[0].ChildNodes.Count; i++) 
                comboBox3.Items.Add(doc.ChildNodes[0].ChildNodes[i].Name);
            return true;
        }
        bool LoadProg (int nom) {
            XmlDocument doc = new XmlDocument();
            try {
                doc.Load(Environment.CurrentDirectory + "\\Save.xml");
            }
            catch {
                return false;
            }
            listBox1.Items.Clear();
            for (int i = 0; i < doc.ChildNodes[0].ChildNodes[nom].ChildNodes.Count; i++) {
                string comand = "";
                string[] buf;
                switch (doc.ChildNodes[0].ChildNodes[nom].ChildNodes[i].Attributes[0].Value) {
                    case "Кнопки_мыши":
                        comand = "Тип команды: Кнопки мыши, команда: ";
                        buf = doc.ChildNodes[0].ChildNodes[nom].ChildNodes[i].Attributes[1].Value.Split('_');
                        comand = comand + buf[0] + " " + buf[1] + " " + buf[2];
                        break;
                    case "Перемешение_мыши":
                        comand = "Тип команды: Перемешение мыши, команда: ";
                        buf = doc.ChildNodes[0].ChildNodes[nom].ChildNodes[i].Attributes[1].Value.Split('_');
                        comand = comand + buf[0] + " " + buf[1];
                        break;
                    case "Клавиатура":
                        comand = "Тип команды: Клавиатура, команда: ";
                        buf = doc.ChildNodes[0].ChildNodes[nom].ChildNodes[i].Attributes[1].Value.Split('_');
                        for (int j = 0; j < buf.Length; j++) 
                            comand = comand + buf[j] + " ";
                        comand = comand.Remove(comand.Length - 1);
                        break;
                    case "Ожидание":
                        comand = "Тип команды: Ожидание, команда: " + doc.ChildNodes[0].ChildNodes[nom].ChildNodes[i].Attributes[1].Value + " сек";
                        break;
                }
                listBox1.Items.Add(comand);
            }
            numericUpDown2.Maximum = listBox1.Items.Count;
            return true;
        }
        void Save () {
            XmlDocument doc = new XmlDocument();
            try {
                doc.Load(Environment.CurrentDirectory + "\\Save.xml");
            }
            catch {
                XmlElement ListXml = doc.CreateElement("Saves");
                doc.AppendChild(ListXml);
            }
            XmlElement progXml = doc.CreateElement(comboBox3.Text);
            doc.ChildNodes[0].AppendChild(progXml);
            for (int i = 0; i < listBox1.Items.Count; i++) {
                XmlElement PoligonXml = doc.CreateElement("Действие_" + i.ToString());
                string[] split = listBox1.Items[i].ToString().Split(' ');
                switch (split[2] + " " + split[3]) {
                    case "Кнопки мыши,":
                        PoligonXml.SetAttribute("Тип_команды", "Кнопки_мыши");
                        PoligonXml.SetAttribute("Команда", split[5] + "_" + split[6] + "_" + split[7]);
                        break;
                    case "Перемешение мыши,":
                        PoligonXml.SetAttribute("Тип_команды", "Перемешение_мыши");
                        PoligonXml.SetAttribute("Команда", split[5] + "_" + split[6]);
                        break;
                    case "Клавиатура, команда:":
                        string word = "";
                        for (int j = 0; j < split.Length - 4; j++)
                            word = word + split[4 + j] + "_";
                        word = word.Remove(word.Length - 1);
                        PoligonXml.SetAttribute("Тип_команды", "Клавиатура");
                        PoligonXml.SetAttribute("Команда", word);
                        break;
                    case "Ожидание, команда:":
                        PoligonXml.SetAttribute("Тип_команды", "Ожидание");
                        PoligonXml.SetAttribute("Команда", split[4]);
                        break;
                }
                progXml.AppendChild(PoligonXml);
            }
            doc.Save(Environment.CurrentDirectory + "\\Save.xml");
        }

        private void button4_Click (object sender, EventArgs e) {
            Save();
            listBox1.Items.Clear();
            comboBox3.SelectedIndex = -1;
            comboBox3.Text = "Новая_программа";
            LoadName();
        }
        private void comboBox3_TextChanged (object sender, EventArgs e) {
            button4.Enabled = !comboBox3.Text.Contains(" ");
            if (comboBox3.Text.Length != 0)
                button4.Enabled = !comboBox3.Text.Contains(" ");
            else
                button4.Enabled = false;
        }

        private void comboBox3_SelectedIndexChanged (object sender, EventArgs e) {
            LoadProg(comboBox3.SelectedIndex);
            size = listBox1.Items.Count;
            if (listBox1.Items.Count > 0) 
                button3.Enabled = button1.Enabled = true;
        }

        private void listBox1_SelectedIndexChanged (object sender, EventArgs e) {
            numericUpDown2.Value = listBox1.SelectedIndex + 1;
        }

        private void numericUpDown2_ValueChanged (object sender, EventArgs e) {
                listBox1.SelectedIndex = (int)numericUpDown2.Value - 1;
        }
    }
}