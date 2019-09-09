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
    partial class Form1  {
        string[] Artist = new string[1];
        public string Folders()  {
            if (folderBrowserDialog1.ShowDialog() == DialogResult.OK) { }
            return folderBrowserDialog1.SelectedPath;
        }

        public void RenameList()  {
            string[] name = Directory.GetFiles(textBox1.Text, "*.mp3");
            name = Rand(name);
            for (int i = 0; i < name.Length; i++)  {
                RenameFile(name[i], i, name.Length);
                progressBar1.Value = (i == (name.Length - 1)) ? 100 : (int)((float)(i + 1F) / (float)(name.Length / 100F));
            }
        }

        string[] Rand(string[] name)  {
            var random = new Random(DateTime.Now.Millisecond);
            return (name.OrderBy(x => random.Next()).ToArray());
        }

        public void RenameFile(string name, int nomer, int Lmas)  {
            string curname = name.Remove(0, name.LastIndexOf('\\') + 1);
            string curpath = name.Remove(name.LastIndexOf('\\') + 1);

            for (int i = 1; i <= 3; i++)
                if (curname[i] == '.')
                    curname = curname.Remove(0, i + 1);

            if (checkBox1.Checked)
                curname = (nomer + 1).ToString() + '.' + curname;

            if (checkBox3.Checked)  {
                int buf = 0;
                for (int i = 0; i < curname.Length; i++)  {
                    if (curname[i] == '-')  {
                        if (curname[i+2] != ' ')
                            buf = i + 2;
                        else
                            buf = i + 3;
                        break;
                    }
                }
                curname = curname.Remove(0, buf);
            }

            var filesmp3 = TagLib.File.Create(name);
            filesmp3.Tag.Title = curname.Remove(curname.Length - 4);
            filesmp3.Tag.Track = (uint)(nomer + 1);
            filesmp3.Tag.TrackCount = (uint)Lmas;
            if (checkBox2.Checked)
                filesmp3.Tag.Performers = Artist;
            filesmp3.Save();
            File.Move(name, (curpath + curname));
        }

        void prov()  {
            bool varprov = false;
            if (Directory.Exists(textBox1.Text))  {
                string[] name = Directory.GetFiles(textBox1.Text, "*.mp3");
                if (name.Length > 1)
                    varprov = true;
            }

            button2.Enabled = varprov;
            checkBox1.Enabled = varprov;
            checkBox3.Enabled = varprov;
            checkBox2.Enabled = varprov;
            if (!varprov)
                checkBox2.Checked = varprov;

            label1.Enabled = checkBox2.Checked;
            textBox2.Enabled = checkBox2.Checked;
            textBox3.Enabled = checkBox2.Checked;
            numericUpDown1.Enabled = checkBox2.Checked;
        }
    }
}