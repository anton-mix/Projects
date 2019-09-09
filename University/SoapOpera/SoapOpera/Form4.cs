using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Forms;
using System.IO;
using System.Windows.Threading;
using System.Drawing;
using System.Runtime.InteropServices;

using Microsoft.Office.Interop.Excel;
using System.Xml;

using VkNet;
using VkNet.Exception;
using VkNet.Enums.Filters;
using VkNet.Infrastructure;
using VkNet.Model;
using VkNet.Model.RequestParams;

namespace SoapOpera {
	public partial class MainWindow : System.Windows.Window {
		DateTime TimeFrom2, TimeTo2;
		private void DeleteComments_Click(object sender, RoutedEventArgs e) {
			if (!ValidForm3())
				return;

			lErrors4.Document.Blocks.Clear();

			var api3 = new VkApi();

			if (!api3.IsAuthorized) {
				try {
					api3.Authorize(new ApiAuthParams {
						ApplicationId = 5539924,
                        Login = ULogin,
                        Password = UserPassword,
                        //AccessToken = "c2d23fe54097f46ea81f8103efd840bc5ca545408c14ee8873c9ee83711d123ca2d5a049f57b774028dfd",
                        Settings = Settings.All,
					});
				}
				catch {
					lErrors4.Document.Blocks.Add(new Paragraph(new Run("Не удалось авторитизироваться в VK. Возможно нет подключения к интеренту.")));
					return;
				}
			}

			try {
				MyFunctions.deleteComments(cbAlbum2.SelectedIndex, AllAlbums, api3, TimeFrom2, TimeTo2);
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Комментарии удалены.")));
			}
			catch (OutOfBorderPhotoException ex) {
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Выход за границу массива фотографий. Критическая ошибка. Exception: " + ex.exception)));
			}
			catch(ErrorFindPhotoException ex) {
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Не удалось получить все фотографии из альбома " + ex.TitleAlbum + ". Возможно отсутствует подключение к интернету. Exception: " + ex.exception)));
			}
			catch(ErrorFindCommentsException ex) {
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Не удалось получить коментарии к фотографии № " + ex.NumberPhoto + " из альбома " + ex.TitleAlbum + ". Возможно отсутствует подключение к интернету. Exception: " + ex.exception)));
			}
			catch (Exception ex){
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Критическая ошибка. " + ex.ToString())));
			}
		}

		private bool ValidForm3() {
			lErrors4.Document.Blocks.Clear();

			if(tbDataFrom2.Text.Length == 0) {
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Не указана дата с которой удалять комментарии")));
				return false;
			}

			if (tbDataTo2.Text.Length == 0) {
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Не указана дата по которую удалять комментарии")));
				return false;
			}

			try {
				string[] splitString = tbDataFrom2.Text.Split(new Char[] { '.' }, StringSplitOptions.RemoveEmptyEntries);
				if (splitString.Length == 3)
					TimeFrom2 = new DateTime(Convert.ToInt32(splitString[2]), Convert.ToInt32(splitString[1]), Convert.ToInt32(splitString[0]));
				else
					throw new Exception();
			}
			catch {
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Не правильно введена дата, с которой начинается поиск. Формат dd.mm.yyyy")));
				return false;
			}

			try {
				string[] splitString = tbDataTo2.Text.Split(new Char[] { '.' }, StringSplitOptions.RemoveEmptyEntries);
				if (splitString.Length == 3)
					TimeTo2 = new DateTime(Convert.ToInt32(splitString[2]), Convert.ToInt32(splitString[1]), Convert.ToInt32(splitString[0]), 23, 59, 59);
				else
					throw new Exception();
			}
			catch {
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Не правильно введена дата, до которой идёт поиск. Формат dd.mm.yyyy")));
				return false;
			}

			if (TimeFrom2 > TimeTo2) {
				lErrors4.Document.Blocks.Add(new Paragraph(new Run("Дата, до которой идёт поиск, должна быть больше даты начала поиска.")));
				return false;
			}

			return true;
		}
	}
}
