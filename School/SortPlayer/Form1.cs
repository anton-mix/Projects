using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;

namespace sortplaer3  {
    public partial class Form1 : Form  {
        public Form1()  {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)  {
            prov();
        }

        private void button1_Click(object sender, EventArgs e)  {
            textBox1.Text = Folders();
            prov();
        }

        private void button2_Click(object sender, EventArgs e)  {
            RenameList();
        }

        private void textBox1_TextChanged(object sender, EventArgs e)  {
            prov();
        }

        private void checkBox2_CheckedChanged(object sender, EventArgs e)  {
            prov();
        }

        private void textBox2_TextChanged(object sender, EventArgs e)  {
            int pr;
            if (Int32.TryParse(textBox2.Text, out pr))  {
                numericUpDown1.Maximum = pr;
                Array.Resize<string>(ref Artist, pr);
            }
            else  {
                Array.Resize<string>(ref Artist, 1);
                numericUpDown1.Maximum = 1;
            }
        }
        private void textBox3_TextChanged(object sender, EventArgs e)  {
            Artist[(int)(numericUpDown1.Value - 1)] = textBox3.Text;
        }

        private void numericUpDown1_ValueChanged(object sender, EventArgs e)  {
            textBox3.Text = Artist[(int)(numericUpDown1.Value - 1)];
        }
    }
}